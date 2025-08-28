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
    visited: bool,

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

    fn newellPlane(face: *const Face, vertices: []const Vertex) [4]f32 {
        var normal: [3]f32 = .{0} ** 3;
        var ctr: [3]f32 = .{0} ** 3;
        var edge = face.edge;
        var n: f32 = 0;
        while (true) {
            const v0 = vertices[edge.tail_vertex];
            const v1 = vertices[edge.next.tail_vertex];
            ctr = add(ctr, v0);
            normal = add(normal, .{
                (v0[1] - v1[1]) * (v0[2] + v1[2]),
                (v0[2] - v1[2]) * (v0[0] + v1[0]),
                (v0[0] - v1[0]) * (v0[1] + v1[1]),
            });
            n += 1;
            edge = edge.next;
            if (edge == face.edge) break;
        }
        const norm = 1.0 / length(normal);
        normal = scale(normal, 1.0 / norm);
        ctr = scale(ctr, 1.0 / n);
        return .{
            normal[0],
            normal[1],
            normal[2],
            -dot(normal, ctr),
        };
    }

    fn addConflict(face: *Face, arena: std.mem.Allocator, vertex: u32) !void {
        const conflict = try arena.create(Conflict);
        conflict.* = .{ .next = face.conflicts, .vertex = vertex };
        face.conflicts = conflict;
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

    // while there are still points outside the polytope, expand it to include them
    while (findConflict(vertices, faces)) |conflict| {
        std.debug.print("{any} {any}\n", .{ conflict.face.plane, vertices[conflict.vertex] });
        conflict.face.visited = true;

        var horizon: std.SegmentedList(*HalfEdge, 4) = .{};
        try buildHorizon(
            arena,
            conflict.face,
            conflict.face.edge,
            vertices[conflict.vertex],
            epsilon,
            &horizon,
        );

        var outside_vertices = try getOutside(arena, faces);
        faces = removeVisited(faces);

        faces = try buildNewFaces(
            arena,
            vertices,
            &horizon,
            faces,
            conflict.vertex,
            epsilon,
            &outside_vertices,
        );
        debugPrint(faces);

        break;
    }

    return faces;
}

fn buildNewFaces(
    arena: std.mem.Allocator,
    vertices: []const Vertex,
    horizon: *std.SegmentedList(*HalfEdge, 4),
    faces: *Face,
    ix_eye: u32,
    epsilon: f32,
    outside_vertices: *std.SegmentedList(u32, 0),
) !*Face {
    var root = faces;
    var new_faces: std.ArrayList(*Face) = try .initCapacity(arena, horizon.count());

    std.debug.print("horizon (len={})\n", .{horizon.count()});
    var it_horizon = horizon.iterator(0);
    while (it_horizon.next()) |h| {
        const edge = h.*;
        const twin = edge.twin;
        std.debug.print("  {}->{}\n", .{ edge.tail_vertex, twin.tail_vertex });

        const new_face = try arena.create(Face);
        new_faces.appendAssumeCapacity(new_face);

        const new_face_edges: [3]*HalfEdge = .{
            try arena.create(HalfEdge),
            try arena.create(HalfEdge),
            try arena.create(HalfEdge),
        };
        // link up to form a circle around tne new face
        for (0..3) |i| {
            new_face_edges[i].next = new_face_edges[(i + 1) % 3];
            new_face_edges[(i + 1) % 3].prev = new_face_edges[i];
            new_face_edges[i].face = new_face;
        }
        // add the twin on the edge touching the horizon
        new_face_edges[0].twin = twin;
        twin.twin = new_face_edges[0];
        new_face_edges[0].tail_vertex = twin.next.tail_vertex;
        new_face_edges[1].tail_vertex = twin.tail_vertex;
        new_face_edges[2].tail_vertex = ix_eye;

        new_face.visited = false;
        new_face.plane = getPlane(
            vertices[new_face_edges[0].tail_vertex],
            vertices[new_face_edges[1].tail_vertex],
            vertices[new_face_edges[2].tail_vertex],
        );
        new_face.edge = new_face_edges[0];

        new_face.conflicts = null;
        var it_outside = outside_vertices.iterator(0);
        while (it_outside.next()) |outside| {
            if (signedDistancePointPlane(vertices[outside.*], new_face.plane) > epsilon) {
                try new_face.addConflict(arena, outside.*);
            }
        }
    }

    // stitch the new faces together
    for (0..new_faces.items.len) |i| {
        const a = new_faces.items[i];
        const b = new_faces.items[(i + 1) % new_faces.items.len];
        a.edge.next.twin = b.edge.prev;
        b.edge.prev.twin = a.edge.next;

        a.prev = null;
        a.next = root;
        root.prev = a;
        root = a;
    }

    return root;
}

// fn doRepairs(...) {

//         merge_loop: while (true) {
//             walk = faces;
//             while (walk) |face| {
//                 var edge = face.edge;
//                 while (true) {
//                     if (!isConvex(vertices, edge.face, edge.twin.face)) {
//                         const new_plane = bestPlane(vertices, edge.face, edge.twin.face);

//                         if (edge.twin.face.prev) |p| {
//                             p.next = edge.twin.face.next;
//                         } else {
//                             faces = edge.twin.face.next.?;
//                         }
//                         if (edge.twin.face.next) |n| {
//                             n.prev = edge.twin.face.prev;
//                         }

//                         face.edge = edge.prev;
//                         face.plane = new_plane;
//                         edge.twin.prev.face = edge.face;
//                         edge.twin.next.face = edge.face;
//                         edge.prev.next = edge.twin.next;
//                         edge.next.prev = edge.twin.prev;
//                         edge.twin.prev.next = edge.next;
//                         edge.twin.next.prev = edge.prev;

//                         topology_loop: while (true) {
//                             edge = face.edge;
//                             while (true) {
//                                 // i'm unsure if we might need to do multiple rounds of repair
//                                 if (edge.twin.face == edge.next.twin.face) {
//                                     face.plane = bestPlane(vertices, face, edge.twin.face);
//                                     if (numSides(edge.twin.face) == 3) {
//                                         // keep the last edge of the twin face
//                                         // drop edge and edge.next and the twin face
//                                         // TODO test correctness
//                                         const last_edge = edge.twin.next;
//                                         last_edge.face = face;
//                                         edge.prev.next = last_edge;
//                                         last_edge.prev = edge.prev;
//                                         edge.next.next.prev = last_edge;
//                                         last_edge.next = edge.next.next;
//                                     } else {
//                                         // keep edge and twin face
//                                         // drop edge.next
//                                         // TODO test correctness
//                                         edge.next.next.prev = edge;
//                                         edge.next = edge.next.next;
//                                         edge.twin.prev.prev = edge.twin;
//                                         edge.twin.prev = edge.twin.prev.prev;
//                                     }

//                                     continue :topology_loop;
//                                 }
//                                 edge = edge.next;
//                                 if (edge == face.edge) break;
//                             }
//                             break;
//                         }

//                         continue :merge_loop;
//                     }

//                     edge = edge.next;
//                     if (edge == face.edge) break;
//                 }
//                 walk = face.next;
//             }

//             break;
//         }

// }

fn getOutside(
    arena: std.mem.Allocator,
    faces: *Face,
) !std.SegmentedList(u32, 0) {
    // this might be N^2 behaviour, it should be valid to check only the outside points
    // of the faces inside the horizon (which might lead to the conflict list being incomplete)
    // however, this is a) simpler and b) seems like it might be more reliable? i dunno
    // maybe come back to this if we need to boost performance
    // simplest pseudo-test for performance would be to just append if visited == true
    var outside_vertices: std.SegmentedList(u32, 0) = .{};
    var walk: ?*Face = faces;
    while (walk) |face| : (walk = face.next) {
        var conflicts: ?*Conflict = face.conflicts;
        while (conflicts) |conflict| : (conflicts = conflict.next) {
            try outside_vertices.append(arena, conflict.vertex);
        }
    }
    return outside_vertices;
}

fn removeVisited(faces: *Face) *Face {
    var root = faces;
    var walk: ?*Face = faces;
    while (walk) |face| {
        const next = face.next;
        if (face.visited) {
            if (face.prev) |p| {
                p.next = next;
            } else {
                root = next.?; // safe since no point can see every face, at least one remains
            }
            if (next) |n| n.prev = face.prev;
        }
        walk = next;
    }
    return root;
}

fn buildHorizon(
    arena: std.mem.Allocator,
    face: *Face,
    starting_edge: *HalfEdge,
    eye: Vertex,
    epsilon: f32,
    horizon: *std.SegmentedList(*HalfEdge, 4),
) !void {
    var edge = starting_edge;
    while (true) {
        const nb = edge.twin.face;
        std.debug.assert(nb != face);
        if (!nb.visited) {
            if (signedDistancePointPlane(eye, nb.plane) < -epsilon) {
                try horizon.append(arena, edge);
            } else {
                nb.visited = true;
                try buildHorizon(arena, nb, edge.twin, eye, epsilon, horizon);
            }
        }
        edge = edge.next;
        if (edge == starting_edge) break;
    }
}

/// walk through the conflict list of every face and find the one farthest away from its face
fn findConflict(vertices: []const Vertex, faces: *Face) ?struct { face: *Face, vertex: u32 } {
    var dist: f32 = -std.math.inf(f32);
    var f: ?*Face = null;
    var ix: u32 = undefined;

    var walk: ?*Face = faces;
    while (walk) |face| : (walk = face.next) {
        var conflicts = face.conflicts;
        while (conflicts) |conflict| : (conflicts = conflict.next) {
            const d = signedDistancePointPlane(vertices[conflict.vertex], face.plane);
            if (d <= dist) continue;
            dist = d;
            f = face;
            ix = conflict.vertex;
        }
    }

    if (f == null) return null;
    return .{ .face = f.?, .vertex = ix };
}

/// generate a large initnial tetrahedron to minimize extra work
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

    face0.visited = false;
    face1.visited = false;
    face2.visited = false;
    face3.visited = false;

    try assertValid(arena, face0, vertices, epsilon);
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

fn dot(a: [3]f32, b: [3]f32) f32 {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
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

    const N = 1;
    for (0..N) |_| {
        const q = qRandom(rand);
        var rotated = try std.ArrayList([3]f32).initCapacity(arena, vertices.len);
        for (vertices) |vertex| rotated.appendAssumeCapacity(qRotate(vertex, q));

        const faces = try quickhull(arena, rotated.items, 1e-6);
        debugPrint(faces);
        try assertValid(arena, faces, rotated.items, 1e-6);
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

fn debugPrint(head: *Face) void {
    std.debug.print("--- begin ---\n", .{});
    var walk: ?*Face = head;
    while (walk) |face| : (walk = face.next) {
        var edge = face.edge;
        std.debug.print("{any}\n  ", .{face.plane});
        while (true) {
            std.debug.print("{}->", .{edge.tail_vertex});
            edge = edge.next;
            if (edge == face.edge) break;
        }
        std.debug.print("{}\n", .{edge.tail_vertex});
    }
    std.debug.print("--- end ---\n", .{});
}

fn assertValid(
    arena: std.mem.Allocator,
    head: *Face,
    vertices: []const Vertex,
    epsilon: f32,
) !void {
    var total_vertex_counts = try arena.alloc(u32, vertices.len);
    for (0..total_vertex_counts.len) |i| total_vertex_counts[i] = 0;

    var local_vertex_counts = try arena.alloc(u32, vertices.len);

    var vertex_pairs: std.AutoHashMapUnmanaged(struct { u32, u32 }, void) = .empty;
    try vertex_pairs.ensureTotalCapacity(arena, @intCast(6 * vertices.len));

    var faces_list: std.AutoHashMapUnmanaged(*Face, void) = .empty;
    try faces_list.ensureTotalCapacity(arena, @intCast(2 * vertices.len));

    // the way this is filled could be more exhaustive
    var faces_dfs: std.AutoHashMapUnmanaged(*Face, void) = .empty;
    try faces_dfs.ensureTotalCapacity(arena, @intCast(2 * vertices.len));

    var n_faces: u32 = 0;
    var n_half_edges: u32 = 0;
    var n_vertices: u32 = 0;

    var walk: ?*Face = head;
    while (walk) |face| {
        var edge = face.edge;
        var n_sides: u32 = 0;

        std.debug.assert(!faces_list.contains(face));
        try faces_list.put(arena, face, {});

        for (0..local_vertex_counts.len) |i| local_vertex_counts[i] = 0;

        while (true) {
            std.debug.assert(edge.next.prev == edge);
            std.debug.assert(edge.prev.next == edge);
            std.debug.assert(edge.twin != edge);
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
            // scale epsilon by the size of the non-normalized dot product
            const u = sub(vertices[edge.next.tail_vertex], vertices[edge.tail_vertex]);
            const v = sub(vertices[edge.prev.tail_vertex], vertices[edge.tail_vertex]);
            std.debug.assert(
                dot(cross(u, v), face.plane[0..3].*) >= -epsilon * length(u) * length(v),
            );

            // the centroid of negihbouring faces must be behind this face
            std.debug.assert(signedDistancePointPlane(
                edge.twin.face.centroid(vertices),
                face.plane,
            ) < epsilon); // since the merge criterion is the opposite of this, i think this holds

            // walk in a cycle around the head vertex
            var walk2 = edge;
            var steps: u32 = 0;
            while (true) {
                walk2 = walk2.twin.next;
                steps += 1;
                std.debug.assert(steps < @as(u32, @intCast(vertices.len * 2)));
                if (walk2 == edge) break;
            }
            std.debug.assert(steps >= 3);

            std.debug.assert(length(sub(
                vertices[edge.next.tail_vertex],
                vertices[edge.tail_vertex],
            )) > epsilon);

            n_sides += 1;
            total_vertex_counts[edge.tail_vertex] += 1;
            local_vertex_counts[edge.tail_vertex] += 1;

            try faces_dfs.put(arena, edge.face, {});
            try faces_dfs.put(arena, edge.twin.face, {});

            std.debug.assert(!vertex_pairs.contains(.{ edge.tail_vertex, edge.next.tail_vertex }));
            try vertex_pairs.put(arena, .{ edge.tail_vertex, edge.next.tail_vertex }, {});

            // all vertices in an edge loop must be unique
            for (local_vertex_counts) |count| std.debug.assert(count < 2);

            edge = edge.next;
            std.debug.assert(n_sides < @as(u32, @intCast(vertices.len * 2)));
            if (edge == face.edge) break;
        }

        // ensure the plane we have isn't corrupted or old
        std.debug.assert(dot(face.newellPlane(vertices)[0..3].*, face.plane[0..3].*) > 0.0);

        n_faces += 1;
        n_half_edges += n_sides;

        std.debug.assert(n_sides >= 3);

        std.debug.assert(n_faces < @as(u32, @intCast(vertices.len * 2)));
        walk = face.next;
    }

    // all vertices (on the hull) must have at least three neighbouring faces
    for (total_vertex_counts) |count| {
        std.debug.assert((count == 0) or (count >= 3));
        if (count > 0) n_vertices += 1;
    }

    // eulers formula
    std.debug.assert(n_half_edges % 2 == 0);
    std.debug.assert(n_vertices + n_faces - n_half_edges / 2 == 2);

    // all faces in the dfs must be in the linked list and vice versa
    std.debug.assert(faces_list.count() == faces_dfs.count());
    var it_face = faces_list.keyIterator();
    while (it_face.next()) |face| std.debug.assert(faces_dfs.contains(face.*));
    it_face = faces_dfs.keyIterator();
    while (it_face.next()) |face| std.debug.assert(faces_list.contains(face.*));
}
