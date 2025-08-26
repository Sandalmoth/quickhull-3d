const std = @import("std");

const log = std.log.scoped(.quickhull_3d);

const Vertex = [3]f32;

const HalfEdge = struct {
    next: *HalfEdge,
    prev: *HalfEdge,
    twin: *HalfEdge,

    tail_vertex: u32,
    face: *Face,
};

const Face = struct {
    next: ?*Face,
    prev: ?*Face,

    plane: [4]f32, // ax + by + cz + d == 0,  a*a + b*b + c*c == 1
    edge: *HalfEdge,
    conflicts: ?*Conflict,

    pub fn centroid(face: *const Face, vertices: []const Vertex) Vertex {
        var edge = face.edge;
        var ctr: [3]f32 = .{0} ** 3;
        var n: f32 = 0;
        while (true) {
            ctr = add(ctr, vertices[edge.tail_vertex]);
            n += 1;
            edge = edge.next;
            if (edge == face.edge) break;
        }
        return scale(ctr, 1.0 / n);
    }

    fn addConflict(face: *Face, arena: std.mem.Allocator, vertex: u32) !void {
        const conflict = try arena.create(Conflict);
        conflict.* = .{ .next = face.conflicts, .vertex = vertex };
        face.conflicts = conflict;
    }

    fn assertValid(head: *Face, vertices: []const Vertex, epsilon: f32) void {
        var walk: ?*Face = head;
        while (walk) |face| {
            var edge = face.edge;
            while (true) {
                std.debug.assert(edge.next.prev == edge);
                std.debug.assert(edge.prev.next == edge);
                std.debug.assert(edge.twin.twin == edge);

                std.debug.assert(edge.face == face);
                std.debug.assert(edge.twin.face != edge.face);

                std.debug.assert(edge.tail_vertex < vertices.len);
                std.debug.assert(edge.twin.tail_vertex == edge.next.tail_vertex);
                std.debug.assert(edge.twin.next.tail_vertex == edge.tail_vertex);

                // all vertices of the plane must be in the plane (within tolerance)
                std.debug.assert(@abs(signedDistancePointPlane(
                    vertices[edge.tail_vertex],
                    face.plane,
                )) < epsilon);

                // face edge loop must be convex
                // not sure what epsilon value to use here though
                // since the given one is interpreted in world space, and this is not
                const u = sub(vertices[edge.next.tail_vertex], vertices[edge.tail_vertex]);
                const v = sub(vertices[edge.prev.tail_vertex], vertices[edge.tail_vertex]);
                std.debug.assert(dot(cross(u, v), face.plane[0..3].*) >= 0); // FIXME epsilon?

                // the centroid of negihbouring faces must be behind this face
                std.debug.assert(signedDistancePointPlane(
                    centroid(edge.twin.face, vertices),
                    face.plane,
                ) < epsilon);

                // walk in a cycle around the tail vertex
                var walk2 = edge;
                while (true) {
                    walk2 = walk2.twin.next;
                    if (walk2 == edge) break;
                }

                edge = edge.next;
                if (edge == face.edge) break;
            }
            walk = face.next;
        }
    }
};

const Conflict = struct {
    next: ?*Conflict,
    vertex: u32,
};

pub fn quickhull(
    arena: std.mem.Allocator,
    vertices: []const Vertex,
    epsilon: f32,
) !*Face {
    var faces = try buildInitialTetrahedron(arena, vertices, epsilon);
    faces = faces;
    return faces;
}

fn buildInitialTetrahedron(
    arena: std.mem.Allocator,
    vertices: []const Vertex,
    epsilon: f32,
) !*Face {
    // a tetrahedron must have at least four faces
    if (vertices.len < 4) return error.Degenerate;

    // first find the axis with the greatest extent and the min/max vertices of that axis
    const extrema: [2]u32 = blk: {
        var extrema: [6]u32 = .{0} ** 6;
        for (vertices[1..], 1..) |vertex, i| {
            if (vertex[0] < vertices[extrema[0]][0]) extrema[0] = @intCast(i);
            if (vertex[1] < vertices[extrema[1]][1]) extrema[1] = @intCast(i);
            if (vertex[2] < vertices[extrema[2]][2]) extrema[2] = @intCast(i);
            if (vertex[0] > vertices[extrema[3]][0]) extrema[3] = @intCast(i);
            if (vertex[1] > vertices[extrema[4]][1]) extrema[4] = @intCast(i);
            if (vertex[2] > vertices[extrema[5]][2]) extrema[5] = @intCast(i);
        }
        var spans: [3]f32 = undefined;
        for (0..3) |i| {
            spans[i] = vertices[extrema[i + 3]][i] - vertices[extrema[i]][i];
        }
        var ix: u32 = 0;
        for (1..3) |i| {
            if (spans[i] > spans[ix]) ix = @intCast(i);
        }
        // all the points are too close together to work at this precision
        if (spans[ix] < epsilon) return error.Degenerate;

        break :blk .{ extrema[ix], extrema[ix + 3] };
    };

    // then find the vertex the greatest distance from the line formed by those vertices
    const furthest: u32 = blk: {
        var dist: f32 = distancePointLine(vertices[0], vertices[extrema[0]], vertices[extrema[1]]);
        var ix: usize = 0;
        for (vertices[1..], 1..) |vertex, i| {
            const d = distancePointLine(vertex, vertices[extrema[0]], vertices[extrema[1]]);
            if (d > dist) {
                ix = @intCast(i);
                dist = d;
            }
        }
        // all the points are too close together to work at this precision
        if (dist < epsilon) return error.Degenerate;
        break :blk @intCast(ix);
    };

    // finally find the vertex farthest away from the plane defined by the previous three
    const final: u32 = blk: {
        const plane = getPlane(vertices[extrema[0]], vertices[extrema[1]], vertices[furthest]);
        var dist: f32 = @abs(signedDistancePointPlane(vertices[0], plane));
        var ix: usize = 0;
        for (vertices[1..], 1..) |vertex, i| {
            const d = @abs(signedDistancePointPlane(vertex, plane));
            if (d > dist) {
                ix = @intCast(i);
                dist = d;
            }
        }
        // again, it must be far enough away to not be considered part of this plane at this eps
        if (dist < epsilon) return error.Degenerate;
        break :blk @intCast(ix);
    };

    // the vertices of the initial tetrahedron
    const tet: [4]u32 = .{
        extrema[0],
        extrema[1],
        furthest,
        final,
    };

    // now construct the half-edge datastructure

    // since there's no constraint on the vertex positions, we have to compute the winding order
    // which should be CCW viewed from the outside, equivalent to the last vertex being behind
    const tri0 = adjustWindingOrder(vertices, .{ tet[1], tet[2], tet[3] }, tet[0]);
    const tri1 = adjustWindingOrder(vertices, .{ tet[0], tet[2], tet[3] }, tet[1]);
    const tri2 = adjustWindingOrder(vertices, .{ tet[0], tet[1], tet[3] }, tet[2]);
    const tri3 = adjustWindingOrder(vertices, .{ tet[0], tet[1], tet[2] }, tet[3]);

    const face0: *Face = try arena.create(Face);
    const face1: *Face = try arena.create(Face);
    const face2: *Face = try arena.create(Face);
    const face3: *Face = try arena.create(Face);

    // regular doubly linked list of faces
    face0.prev = null;
    face0.next = face1;
    face1.prev = face0;
    face1.next = face2;
    face2.prev = face1;
    face2.next = face3;
    face3.prev = face2;
    face3.next = null;

    face0.plane = getPlane(vertices[tri0[0]], vertices[tri0[1]], vertices[tri0[2]]);
    face1.plane = getPlane(vertices[tri1[0]], vertices[tri1[1]], vertices[tri1[2]]);
    face2.plane = getPlane(vertices[tri2[0]], vertices[tri2[1]], vertices[tri2[2]]);
    face3.plane = getPlane(vertices[tri3[0]], vertices[tri3[1]], vertices[tri3[2]]);

    const face0_edges: [3]*HalfEdge = .{
        try arena.create(HalfEdge),
        try arena.create(HalfEdge),
        try arena.create(HalfEdge),
    };
    const face1_edges: [3]*HalfEdge = .{
        try arena.create(HalfEdge),
        try arena.create(HalfEdge),
        try arena.create(HalfEdge),
    };
    const face2_edges: [3]*HalfEdge = .{
        try arena.create(HalfEdge),
        try arena.create(HalfEdge),
        try arena.create(HalfEdge),
    };
    const face3_edges: [3]*HalfEdge = .{
        try arena.create(HalfEdge),
        try arena.create(HalfEdge),
        try arena.create(HalfEdge),
    };

    face0.edge = face0_edges[0];
    face1.edge = face1_edges[0];
    face2.edge = face2_edges[0];
    face3.edge = face3_edges[0];

    // edges form a circular linked list around each face
    for (0..3) |i| {
        face0_edges[i].next = face0_edges[(i + 1) % 3];
        face0_edges[(i + 1) % 3].prev = face0_edges[i];
        face0_edges[i].tail_vertex = tri0[i];
        face0_edges[i].face = face0;

        face1_edges[i].next = face1_edges[(i + 1) % 3];
        face1_edges[(i + 1) % 3].prev = face1_edges[i];
        face1_edges[i].tail_vertex = tri1[i];
        face1_edges[i].face = face1;

        face2_edges[i].next = face2_edges[(i + 1) % 3];
        face2_edges[(i + 1) % 3].prev = face2_edges[i];
        face2_edges[i].tail_vertex = tri2[i];
        face2_edges[i].face = face2;

        face3_edges[i].next = face3_edges[(i + 1) % 3];
        face3_edges[(i + 1) % 3].prev = face3_edges[i];
        face3_edges[i].tail_vertex = tri3[i];
        face3_edges[i].face = face3;
    }

    // for each edge find its twin, i.e. the other edge that connects to the same vertices
    // this could be simplified by writing both sides of the twin at the same time
    // could it be simplified further by analyzing the vertex lists to link the right sides?
    for (0..3) |i| {
        for (0..3) |j| {
            if ((face0_edges[i].tail_vertex == face1_edges[j].next.tail_vertex) and
                (face0_edges[i].next.tail_vertex == face1_edges[j].tail_vertex))
                face0_edges[i].twin = face1_edges[j];
            if ((face0_edges[i].tail_vertex == face2_edges[j].next.tail_vertex) and
                (face0_edges[i].next.tail_vertex == face2_edges[j].tail_vertex))
                face0_edges[i].twin = face2_edges[j];
            if ((face0_edges[i].tail_vertex == face3_edges[j].next.tail_vertex) and
                (face0_edges[i].next.tail_vertex == face3_edges[j].tail_vertex))
                face0_edges[i].twin = face3_edges[j];

            if ((face1_edges[i].tail_vertex == face0_edges[j].next.tail_vertex) and
                (face1_edges[i].next.tail_vertex == face0_edges[j].tail_vertex))
                face1_edges[i].twin = face0_edges[j];
            if ((face1_edges[i].tail_vertex == face2_edges[j].next.tail_vertex) and
                (face1_edges[i].next.tail_vertex == face2_edges[j].tail_vertex))
                face1_edges[i].twin = face2_edges[j];
            if ((face1_edges[i].tail_vertex == face3_edges[j].next.tail_vertex) and
                (face1_edges[i].next.tail_vertex == face3_edges[j].tail_vertex))
                face1_edges[i].twin = face3_edges[j];

            if ((face2_edges[i].tail_vertex == face0_edges[j].next.tail_vertex) and
                (face2_edges[i].next.tail_vertex == face0_edges[j].tail_vertex))
                face2_edges[i].twin = face0_edges[j];
            if ((face2_edges[i].tail_vertex == face1_edges[j].next.tail_vertex) and
                (face2_edges[i].next.tail_vertex == face1_edges[j].tail_vertex))
                face2_edges[i].twin = face1_edges[j];
            if ((face2_edges[i].tail_vertex == face3_edges[j].next.tail_vertex) and
                (face2_edges[i].next.tail_vertex == face3_edges[j].tail_vertex))
                face2_edges[i].twin = face3_edges[j];

            if ((face3_edges[i].tail_vertex == face0_edges[j].next.tail_vertex) and
                (face3_edges[i].next.tail_vertex == face0_edges[j].tail_vertex))
                face3_edges[i].twin = face0_edges[j];
            if ((face3_edges[i].tail_vertex == face1_edges[j].next.tail_vertex) and
                (face3_edges[i].next.tail_vertex == face1_edges[j].tail_vertex))
                face3_edges[i].twin = face1_edges[j];
            if ((face3_edges[i].tail_vertex == face2_edges[j].next.tail_vertex) and
                (face3_edges[i].next.tail_vertex == face2_edges[j].tail_vertex))
                face3_edges[i].twin = face2_edges[j];
        }
    }

    face0.conflicts = null;
    face1.conflicts = null;
    face2.conflicts = null;
    face3.conflicts = null;

    for (vertices, 0..) |vertex, j| {
        const i: u32 = @intCast(j);
        if (signedDistancePointPlane(vertex, face0.plane) > epsilon) try face0.addConflict(arena, i);
        if (signedDistancePointPlane(vertex, face1.plane) > epsilon) try face1.addConflict(arena, i);
        if (signedDistancePointPlane(vertex, face2.plane) > epsilon) try face2.addConflict(arena, i);
        if (signedDistancePointPlane(vertex, face3.plane) > epsilon) try face3.addConflict(arena, i);
    }

    face0.assertValid(vertices, epsilon);
    return face0;
}

fn adjustWindingOrder(vertices: []const Vertex, tri: [3]u32, other: u32) [3]u32 {
    const plane = getPlane(vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]);
    if (signedDistancePointPlane(vertices[other], plane) < 0) {
        return tri;
    } else {
        return .{ tri[2], tri[1], tri[0] };
    }
}

fn distancePointLine(p: [3]f32, a: [3]f32, b: [3]f32) f32 {
    const ab = sub(b, a);
    const ap = sub(p, a);
    const c = cross(ap, ab);
    const area = length(c);
    const base = length(ab);
    return if (base > 0) area / base else base;
}

fn getPlane(a: [3]f32, b: [3]f32, c: [3]f32) [4]f32 {
    const u = sub(b, a);
    const v = sub(c, a);
    var n = cross(u, v);
    const norm = 1.0 / length(n);
    for (0..3) |i| n[i] *= norm;
    return .{ n[0], n[1], n[2], -(n[0] * a[0] + n[1] * a[1] + n[2] * a[2]) };
}

fn signedDistancePointPlane(point: [3]f32, plane: [4]f32) f32 {
    return plane[0] * point[0] + plane[1] * point[1] + plane[2] * point[2] + plane[3];
}

fn sub(a: [3]f32, b: [3]f32) [3]f32 {
    return .{
        a[0] - b[0],
        a[1] - b[1],
        a[2] - b[2],
    };
}

fn add(a: [3]f32, b: [3]f32) [3]f32 {
    return .{
        a[0] + b[0],
        a[1] + b[1],
        a[2] + b[2],
    };
}

fn scale(a: [3]f32, x: f32) [3]f32 {
    return .{
        a[0] * x,
        a[1] * x,
        a[2] * x,
    };
}

fn cross(a: [3]f32, b: [3]f32) [3]f32 {
    return .{
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    };
}

fn length(a: [3]f32) f32 {
    return @sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

test "rotated dodecahedron" {
    var arena_impl = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena_impl.deinit();
    const arena = arena_impl.allocator();

    const vertices = [_][3]f32{
        .{ -1, 1, -1 },
        .{ -1, 1, 1 },
        .{ 1, -1, 1 },
        .{ 1, -1, -1 },
        .{ -0.6180339887498948, 0, -1.618033988749895 },
        .{ 0, 1.618033988749895, -0.6180339887498948 },
        .{ 0, -1.618033988749895, -0.6180339887498948 },
        .{ -0.6180339887498948, 0, 1.618033988749895 },
        .{ -1.618033988749895, -0.6180339887498948, 0 },
        .{ 1.618033988749895, 0.6180339887498948, 0 },
        .{ -1, -1, -1 },
        .{ -1, -1, 1 },
        .{ 1.618033988749895, -0.6180339887498948, 0 },
        .{ -1.618033988749895, 0.6180339887498948, 0 },
        .{ 0, 1.618033988749895, 0.6180339887498948 },
        .{ 0, -1.618033988749895, 0.6180339887498948 },
        .{ 0.6180339887498948, 0, -1.618033988749895 },
        .{ 1, 1, -1 },
        .{ 1, 1, 1 },
        .{ 0.6180339887498948, 0, 1.618033988749895 },
    };

    var rng = std.Random.DefaultPrng.init(1337);
    const rand = rng.random();

    const N = 1000;
    for (0..N) |_| {
        const q = qRandom(rand);
        var rotated = try std.ArrayList([3]f32).initCapacity(arena, vertices.len);
        for (vertices) |vertex| rotated.appendAssumeCapacity(qRotate(vertex, q));

        const faces = try quickhull(arena, rotated.items, 1e-5);
        _ = faces;
    }
}

fn qRotate(p: [3]f32, q: [4]f32) [3]f32 {
    var a = cross(.{ q[0], q[1], q[2] }, p);
    for (0..3) |i| a[i] *= 2;
    const b = cross(.{ q[0], q[1], q[2] }, a);
    var p_rot: [3]f32 = undefined;
    for (0..3) |i| p_rot[i] = p[i] + q[3] * a[i] + b[i];
    return p_rot;
}

fn qRandom(rand: std.Random) [4]f32 {
    const x = rand.floatNorm(f32);
    const y = rand.floatNorm(f32);
    const z = rand.floatNorm(f32);
    const w = rand.floatNorm(f32);
    const norm = @sqrt(x * x + y * y + z * z + w * w);
    if (norm < 1e-12) return qRandom(rand);
    const inorm = 1.0 / norm;
    return .{
        inorm * x,
        inorm * y,
        inorm * z,
        inorm * w,
    };
}
fn dot(a: [3]f32, b: [3]f32) f32 {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
