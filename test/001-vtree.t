#!/usr/bin/env escript
%% -*- erlang -*-

% Licensed under the Apache License, Version 2.0 (the "License"); you may not
% use this file except in compliance with the License. You may obtain a copy of
% the License at
%
%   http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
% WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
% License for the specific language governing permissions and limitations under
% the License.

-define(FILENAME, "/tmp/vtree.bin").

main(_) ->
    code:add_pathz(filename:dirname(escript:script_name())),
    gc_test_util:init_code_path(),
    etap:plan(138),
    case (catch test()) of
        ok ->
            etap:end_tests();
        Other ->
            etap:diag(io_lib:format("Test died abnormally: ~p", [Other])),
            etap:bail(Other)
    end,
    ok.

test() ->
    test_intersect(),
    test_disjoint(),
    test_lookup(),
    test_multilookup(),
    test_knn(),
    test_split_bbox_if_flipped(),
    test_area(),
    test_merge_mbr(),
    test_find_area_min_nth(),
    test_partition_node(),
    test_calc_mbr(),
    test_calc_nodes_mbr(),
    test_best_split(),
    test_minimal_overlap(),
    test_minimal_coverage(),
    test_calc_overlap(),
    test_insert(),
    test_delete(),
    test_delete_same_id(),
    test_split_node(),
    test_count_total(),

    %etap:end_tests().
    ok.

-record(node, {
    % type = inner | leaf
    type=leaf}).

test_insert() ->
    %etap:plan(6),

    {ok, Fd} = case couch_file:open(?FILENAME, [create, overwrite]) of
    {ok, Fd2} ->
        {ok, Fd2};
    {error, Reason} ->
        io:format("ERROR (~s): Couldn't open file (~s) for tree storage~n",
                  [Reason, ?FILENAME])
    end,

    Node1 = {{10,5,13,15}, #node{type=leaf},
        {linestring, [[10,5],[13,15]]}, <<"Node1">>},
    Node2 = {{-18,-3,-10,-1}, #node{type=leaf},
        {linestring, [[-18,-3],[-10,-1]]}, <<"Node2">>},
    Node3 = {{-21,2,-10,14}, #node{type=leaf},
        {linestring, [[-21,2],[-10,14]]}, <<"Node3">>},
    Node4 = {{5,-32,19,-25}, #node{type=leaf},
        {linestring, [[5,-32],[19,-25]]}, <<"Node4">>},
    Node5 = {{-5,-16,4,19}, #node{type=leaf},
        {linestring, [[-5,-16],[4,19]]}, <<"Node5">>},
    % od = on disk
    Node1od = {{10,5,13,15}, #node{type=leaf},
        {<<"Node1">>, {{linestring, [[10,5],[13,15]]}, <<"Node1">>}}},
    Node2od = {{-18,-3,-10,-1}, #node{type=leaf},
        {<<"Node2">>, {{linestring, [[-18,-3],[-10,-1]]}, <<"Node2">>}}},
    Node3od = {{-21,2,-10,14}, #node{type=leaf},
        {<<"Node3">>, {{linestring, [[-21,2],[-10,14]]}, <<"Node3">>}}},
    Node4od = {{5,-32,19,-25}, #node{type=leaf},
        {<<"Node4">>, {{linestring, [[5,-32],[19,-25]]}, <<"Node4">>}}},
    Node5od = {{-5,-16,4,19}, #node{type=leaf},
        {<<"Node5">>, {{linestring, [[-5,-16],[4,19]]}, <<"Node5">>}}},
    Mbr1 = {10,5,13,15},
    Mbr1_2 = {-18,-3,13,15},
    Mbr1_2_3 = {-21,-3,13,15},
    Mbr1_2_3_4 = {-21,-32,19,15},
    Mbr1_2_3_4_5 = {-21,-32,19,19},
    Mbr1_4_5 = {-5,-32,19,19},
    Mbr2_3 = {-21,-3,-10,14},
    Tree1Node1 = {Mbr1, #node{type=leaf}, [Node1od]},
    Tree1Node1_2 = {Mbr1_2, #node{type=leaf}, [Node1od, Node2od]},
    Tree1Node1_2_3 = {Mbr1_2_3, #node{type=leaf}, [Node1od, Node2od, Node3od]},
    Tree1Node1_2_3_4 = {Mbr1_2_3_4, #node{type=leaf},
                        [Node1od, Node2od, Node3od, Node4od]},
    Tree1Node1_2_3_4_5 = {Mbr1_2_3_4_5, #node{type=inner},
                          [{ok, {Mbr2_3, #node{type=leaf}, [Node2od, Node3od]}},
                           {ok, {Mbr1_4_5, #node{type=leaf},
                                 [Node1od, Node4od, Node5od]}}]},

    etap:is(vtree:insert(Fd, nil, <<"Node1">>, Node1), {ok, Mbr1, 0, 1},
            "Insert a node into an empty tree (write to disk)"),
    etap:is(get_node(Fd, 0), {ok, Tree1Node1},
            "Insert a node into an empty tree" ++
            " (check if it was written correctly)"),
    {ok, Mbr1_2, Pos2, 1} = vtree:insert(Fd, 0, <<"Node2">>, Node2),
    etap:is(get_node(Fd, Pos2), {ok, Tree1Node1_2},
            "Insert a node into a not yet full leaf node (root node) (a)"),
    {ok, Mbr1_2_3, Pos3, 1} = vtree:insert(Fd, Pos2, <<"Node3">>, Node3),
    etap:is(get_node(Fd, Pos3), {ok, Tree1Node1_2_3},
            "Insert a node into a not yet full leaf node (root node) (b)"),
    {ok, Mbr1_2_3_4, Pos4, 1} = vtree:insert(Fd, Pos3, <<"Node4">>, Node4),
    etap:is(get_node(Fd, Pos4), {ok, Tree1Node1_2_3_4},
            "Insert a nodes into a then to be full leaf node (root node)"),
    {ok, Mbr1_2_3_4_5, Pos5, 2} = vtree:insert(Fd, Pos4, <<"Node5">>, Node5),
    {ok, {Mbr1_2_3_4_5, #node{type=inner}, [Pos5_1, Pos5_2]}} =
                get_node(Fd, Pos5),
    etap:is({ok, {Mbr1_2_3_4_5, #node{type=inner},
                  [get_node(Fd, Pos5_1), get_node(Fd, Pos5_2)]}},
                  {ok, Tree1Node1_2_3_4_5},
            "Insert a nodes into a full leaf node (root node)"),
    ok.


test_intersect() ->
    %etap:plan(17),
    Mbr1_2 = {-18,-3,13,15},
    etap:is(vtree:intersect(Mbr1_2, {-20, -11, 0, 0}), true,
            "MBR intersectton (S and W edge)"),
    etap:is(vtree:intersect(Mbr1_2, {-21, 4, -2, 11}), true,
            "MBR intersecttion (W edge)"),
    etap:is(vtree:intersect(Mbr1_2, {-21, 4, -2, 17}), true,
            "MBR intersection (W and N edge)"),
    etap:is(vtree:intersect(Mbr1_2, {-13, 4, -2, 17}), true,
            "MBR intersection (N edge)"),
    etap:is(vtree:intersect(Mbr1_2, {-13, 4, 16, 17}), true,
            "MBR intersection (N and E edge)"),
    etap:is(vtree:intersect(Mbr1_2, {5, -1, 16, 10}), true,
            "MBR intersection (E edge)"),
    etap:is(vtree:intersect(Mbr1_2, {5, -9, 16, 10}), true,
            "MBR intersection (E and S edge)"),
    etap:is(vtree:intersect(Mbr1_2, {5, -9, 11, 10}), true,
            "MBR intersection (S edge)"),
    etap:is(vtree:intersect(Mbr1_2, {-27, -9, 11, 10}), true,
            "MBR intersection (S and W edge)"),
    etap:is(vtree:intersect(Mbr1_2, {-27, -9, 18, 10}), true,
            "MBR intersection (W and E edge (bottom))"),
    etap:is(vtree:intersect(Mbr1_2, {-27, 2, 18, 10}), true,
            "MBR intersection (W and E edge (middle))"),
    etap:is(vtree:intersect(Mbr1_2, {-10, -9, 18, 10}), true,
            "MBR intersection (W and E edge (top))"),
    etap:is(vtree:intersect(Mbr1_2, {-25, -4, 2, 12}), true,
            "MBR intersection (N and S edge (left))"),
    etap:is(vtree:intersect(Mbr1_2, {-15, -4, 2, 12}), true,
            "MBR intersection (N and S edge (middle))"),
    etap:is(vtree:intersect(Mbr1_2, {-15, -4, 2, 22}), true,
            "MBR intersection (N and S edge (right))"),
    etap:is(vtree:intersect(Mbr1_2, {-14, -1, 10, 5}), false,
            "One MBR within the other"),
    etap:is(vtree:intersect(Mbr1_2, Mbr1_2), true,
            "MBR is within itself"),
    ok.

test_disjoint() ->
    %etap:plan(2),
    Mbr1_2 = {-18,-3,13,15},
    etap:is(vtree:disjoint(Mbr1_2, {27, 20, 38, 40}), true,
            "MBRs are disjoint"),
    etap:is(vtree:disjoint(Mbr1_2, {-27, 2, 18, 10}), false,
            "MBRs are not disjoint").


test_lookup() ->
    %etap:plan(8),

    {ok, Fd} = case couch_file:open(?FILENAME, [create, overwrite]) of
    {ok, Fd2} ->
        {ok, Fd2};
    {error, Reason} ->
        io:format("ERROR (~s): Couldn't open file (~s) for tree storage~n",
                  [Reason, ?FILENAME])
    end,
    Node1 = {{10,5,13,15}, #node{type=leaf},
        {linestring, [[10,5],[13,15]]}, <<"Node1">>},
    Node2 = {{-18,-3,-10,-1}, #node{type=leaf},
        {linestring, [[-18,-3],[-10,-1]]}, <<"Node2">>},
    Node3 = {{-21,2,-10,14}, #node{type=leaf},
        {linestring, [[-21,2],[-10,14]]}, <<"Node3">>},
    Node4 = {{5,-32,19,-25}, #node{type=leaf},
        {linestring, [[5,-32],[19,-25]]}, <<"Node4">>},
    Node5 = {{-5,-16,4,19}, #node{type=leaf},
        {linestring, [[-5,-16],[4,19]]}, <<"Node5">>},
    {Mbr1, _, Geom1, Id1} = Node1,
    {Mbr2, _, Geom2, Id2} = Node2,
    {_, _, _Geom3, Id3} = Node3,
    {_, _, _Geom4, Id4} = Node4,
    {Mbr5, _, Geom5, Id5} = Node5,
    Mbr1_2 = {-18,-3,13,15},
    Mbr1_2_3 = {-21,-3,13,15},
    Mbr1_2_3_4 = {-21,-32,19,15},
    Mbr1_2_3_4_5 = {-21,-32,19,19},
    Bbox1 = {-20, -10, 30, 21},
    Bbox2 = {-20, -10, 0, 0},
    Bbox3 = {100, 200, 300, 400},
    Bbox4 = {0, 0, 20, 15},

    {ok, Mbr1, 0, 1} = vtree:insert(Fd, nil, Id1, Node1),
    {ok, Mbr1_2, Pos2, 1} = vtree:insert(Fd, 0, Id2, Node2),
    {ok, Mbr1_2_3, Pos3, 1} = vtree:insert(Fd, Pos2, Id3, Node3),
    {ok, Mbr1_2_3_4, Pos4, 1} = vtree:insert(Fd, Pos3, Id4, Node4),
    {ok, Mbr1_2_3_4_5, Pos5, 2} = vtree:insert(Fd, Pos4, Id5, Node5),

    {ok, Lookup1} = gc_test_util:lookup(Fd, Pos2, Bbox1),
    etap:is(Lookup1, [{Mbr2, Id2, Geom2, Id2},{Mbr1, Id1, Geom1, Id1}],
            "Find all nodes in tree (tree height=1)"),
    {ok, Lookup2} = gc_test_util:lookup(Fd, Pos2, Bbox2),
    etap:is(Lookup2, [{Mbr2, Id2, Geom2, Id2}],
            "Find some nodes in tree (tree height=1)"),
    {ok, Lookup3} = gc_test_util:lookup(Fd, Pos2, Bbox3),
    etap:is(Lookup3, [], "Query window outside of all nodes (tree height=1)"),
    {ok, Lookup4} = gc_test_util:lookup(Fd, Pos5, Bbox2),
    etap:is(Lookup4, [{Mbr5, Id5, Geom5, Id5}, {Mbr2, Id2, Geom2, Id2}],
            "Find some nodes in tree (tree height=2) (a)"),
    {ok, Lookup5} = gc_test_util:lookup(Fd, Pos5, Bbox4),
    etap:is(Lookup5, [{Mbr5, Id5, Geom5, Id5}, {Mbr1, Id1, Geom1, Id1}],
            "Find some nodes in tree (tree height=2) (b)"),
    {ok, Lookup6} = gc_test_util:lookup(Fd, Pos5, Bbox3),
    etap:is(Lookup6, [], "Query window outside of all nodes (tree height=2)"),

    {ok, Lookup7} = gc_test_util:lookup(Fd, Pos5, [Bbox2, Bbox4]),
    etap:is(length(Lookup7), 3, "Query with multiple windows (2 windows)"),
    {ok, Lookup8} = gc_test_util:lookup(Fd, Pos5, [Bbox2, Bbox4, {-20,1,-9,15}]),
    etap:is(length(Lookup8), 4, "Query with multiple windows (3 windows)"),
    ok.


test_multilookup() ->
    %etap:plan(5),

    {ok, Fd} = case couch_file:open(?FILENAME, [create, overwrite]) of
    {ok, Fd2} ->
        {ok, Fd2};
    {error, Reason} ->
        io:format("ERROR (~s): Couldn't open file (~s) for tree storage~n",
                  [Reason, ?FILENAME])
    end,

    Node1 = {{10,5,13,15}, #node{type=leaf},
        {linestring, [[10,5],[13,15]]}, <<"Node1">>},
    Node2 = {{-18,-3,-10,-1}, #node{type=leaf},
        {linestring, [[-18,-3],[-10,-1]]}, <<"Node2">>},
    Node3 = {{-21,2,-10,14}, #node{type=leaf},
        {linestring, [[-21,2],[-10,14]]}, <<"Node3">>},
    Node4 = {{5,-32,19,-25}, #node{type=leaf},
        {linestring, [[5,-32],[19,-25]]}, <<"Node4">>},
    Node5 = {{-5,-16,4,19}, #node{type=leaf},
        {linestring, [[-5,-16],[4,19]]}, <<"Node5">>},
    {Mbr1, _, Geom1, Id1} = Node1,
    {Mbr2, _, Geom2, Id2} = Node2,
    {_, _, _, Id3} = Node3,
    {_, _, _, Id4} = Node4,
    {_, _, _, Id5} = Node5,
    Mbr1_2 = {-18,-3,13,15},
    Mbr1_2_3 = {-21,-3,13,15},
    Mbr1_2_3_4 = {-21,-32,19,15},
    Mbr1_2_3_4_5 = {-21,-32,19,19},
    Bbox1 = {-20, -10, 30, 21},
    Bbox1a = {-20, -10, 5, 10},
    Bbox1b = {5, 10, 30, 21},
    Bbox2a = {-20, -17, -17, 0},
    Bbox2b = {-15, -5, -8, -2},
    Bbox2c = {200, 100, 190, 90},
    Bbox3 = {100, 200, 300, 400},

    {ok, Mbr1, 0, 1} = vtree:insert(Fd, nil, Id1, Node1),
    {ok, Mbr1_2, Pos2, 1} = vtree:insert(Fd, 0, Id2, Node2),
    {ok, Mbr1_2_3, Pos3, 1} = vtree:insert(Fd, Pos2, Id3, Node3),
    {ok, Mbr1_2_3_4, Pos4, 1} = vtree:insert(Fd, Pos3, Id4, Node4),
    {ok, Mbr1_2_3_4_5, Pos5, 2} = vtree:insert(Fd, Pos4, Id5, Node5),

    {ok, Lookup1} = gc_test_util:lookup(Fd, Pos2, [Bbox1]),
    etap:is(Lookup1, [{Mbr2, Id2, Geom2, Id2}, {Mbr1, Id1, Geom1, Id1}],
            "Find all nodes in tree (tree height=1) with one bbox"),
    {ok, Lookup2} = gc_test_util:lookup(Fd, Pos2, [Bbox1a, Bbox1b]),
    etap:is(Lookup2, [{Mbr2, Id2, Geom2, Id2}, {Mbr1, Id1, Geom1, Id1}],
            "Find all nodes in tree (tree height=1) with two bboxes"),
    {ok, Lookup3} = gc_test_util:lookup(Fd, Pos5, [Bbox2a, Bbox2b]),
    etap:is(Lookup3, [{Mbr2, Id2, Geom2, Id2}],
            "Find some nodes in tree (tree height=2) 2 bboxes"),
    {ok, Lookup4} = gc_test_util:lookup(Fd, Pos5, [Bbox2a, Bbox2b, Bbox2c]),
    etap:is(Lookup4, [{Mbr2, Id2, Geom2, Id2}],
            "Find some nodes in tree (tree height=2) 3 bboxes (one outside "
            "of all nodes"),
    {ok, Lookup5} = gc_test_util:lookup(Fd, Pos5, [Bbox2c, Bbox3]),
    etap:is(Lookup5, [],
            "Don't find any nodes in tree (tree height=2) 2 bboxes (both "
            "outside of all nodes"),
    ok.

test_knn() ->
    %etap:plan(5),

    {ok, Fd} = case couch_file:open(?FILENAME, [create, overwrite]) of
    {ok, Fd2} ->
        {ok, Fd2};
    {error, Reason} ->
        io:format("ERROR (~s): Couldn't open file (~s) for tree storage~n",
                  [Reason, ?FILENAME])
    end,
    Node1 = {{1,1,2,2}, #node{type=leaf},
        {linestring, [[1,1],[2,2]]}, <<"Node1">>},
    Node2 = {{3,1,4,2}, #node{type=leaf},
        {linestring, [[3,2],[4,1]]}, <<"Node2">>},
    Node3 = {{1,4,2,5}, #node{type=leaf},
        {linestring, [[1,5],[2,4]]}, <<"Node3">>},
    Node4 = {{4,5,5,6}, #node{type=leaf},
        {linestring, [[4,5],[5,6]]}, <<"Node4">>},
    Node5 = {{6,3,7,4}, #node{type=leaf},
        {linestring, [[6,3],[7,4]]}, <<"Node5">>},
    Node6 = {{8.5,1,9.5,2}, #node{type=leaf},
        {linestring, [[8.5,2],[9.5,1]]}, <<"Node6">>},
    {Mbr1, _, Geom1, Id1} = Node1,
    {_, _, _, Id2} = Node2,
    {_, _, _, Id3} = Node3,
    {_, _, _, Id4} = Node4,
    {_, _, _, Id5} = Node5,
    {_, _, _, Id6} = Node6,
    Mbr1_2 = {1,1,4,2},
    Mbr1_2_3 = {1,1,4,5},
    Mbr1_2_3_4 = {1,1,5,6},
    Mbr1_2_3_4_5 = {1,1,7,6},
    Mbr1_2_3_4_5_6 = {1,1,9.5,6},

    {ok, Mbr1, 0, 1} = vtree:insert(Fd, nil, Id1, Node1),
    {ok, Mbr1_2, Pos2, 1} = vtree:insert(Fd, 0, Id2, Node2),
    {ok, Mbr1_2_3, Pos3, 1} = vtree:insert(Fd, Pos2, Id3, Node3),
    {ok, Mbr1_2_3_4, Pos4, 1} = vtree:insert(Fd, Pos3, Id4, Node4),
    {ok, Mbr1_2_3_4_5, Pos5, 2} = vtree:insert(Fd, Pos4, Id5, Node5),
    {ok, Mbr1_2_3_4_5_6, Pos6, 2} = vtree:insert(Fd, Pos5, Id6, Node6),

    {ok, Result1} = gc_test_util:knn(Fd, Pos2, 1, {-45, -45}, nil, nil),
    etap:is(Result1, [{Mbr1, Id1, Geom1, Id1}],
            "Get single node with Id, Mbr, Geom and Value"),

    {ok, Result2} = gc_test_util:knnIds(Fd, Pos6, 1000, {0, 0}, nil, nil),
    etap:is(Result2, [Id6, Id5, Id4, Id3, Id2, Id1],
            "Get all nodes, sorted by distance desc."),

    {ok, Result3} = gc_test_util:knnIds(Fd, Pos6, 2, {2.5, 0}, nil, nil),
    etap:is(Result3, [Id2, Id1],
            "Get 2 nodes with same distance"),

    Bounds = {0,0,10,10},
    {ok, Result4} = gc_test_util:knnIds(Fd, Pos6, 1, {0,1}, Bounds, nil),
    etap:is(Result4, [Id6],
            "Bounds are used for distance calculation"),

    {ok, Result5} = gc_test_util:knnIds(Fd, Pos6, 1000, {-45, -45}, nil, true),
    etap:is(Result5, [Id4, Id5, Id6, Id3, Id2, Id1],
            "Get all nodes, sorted by spherical distance desc."),
    ok.


test_split_bbox_if_flipped() ->
    %etap:plan(17),

    etap:is(vtree:split_bbox_if_flipped({160,40,-120,60}, {-180,-90,180,90}),
            [{-180,40,-120,60}, {160,40,180,60}],
            "Bbox over the date line and north (flipped in x direction)"),
    etap:is(vtree:split_bbox_if_flipped({160,-60,-120,-40}, {-180,-90,180,90}),
            [{-180,-60,-120,-40}, {160,-60,180,-40}],
            "Bbox over the date line and south (flipped in x direction)"),
    etap:is(vtree:split_bbox_if_flipped({160,-40,-120,60}, {-180,-90,180,90}),
            [{-180,-40,-120,60}, {160,-40,180,60}],
            "Bbox over the date line and over equator (flipped in x "
            "direction)"),

    etap:is(vtree:split_bbox_if_flipped({20,70,30,-60}, {-180,-90,180,90}),
            [{20,-90,30,-60}, {20,70,30,90}],
            "Bbox over the pole and east (flipped in y direction)"),
    etap:is(vtree:split_bbox_if_flipped({-30,70,-20,-60}, {-180,-90,180,90}),
            [{-30,-90,-20,-60}, {-30,70,-20,90}],
            "Bbox over the pole and west (flipped in y direction)"),
    etap:is(vtree:split_bbox_if_flipped({-30,70,20,-60}, {-180,-90,180,90}),
            [{-30,-90,20,-60}, {-30,70,20,90}],
            "Bbox over the pole and Greenwich (flipped in y direction)"),

    etap:is(vtree:split_bbox_if_flipped({150,70,30,-60}, {-180,-90,180,90}),
            [{-180,-90,30,-60}, {-180,70,30,90},
             {150,-90,180,-60}, {150,70,180,90}],
            "Bbox over the pole and date line (flipped in x and y direction)"),
    etap:is(vtree:split_bbox_if_flipped({-30,-60,-20,70}, {-180,-90,180,90}),
            [{-30,-60,-20,70}],
            "Not flipped Bbox"),

    etap:is(vtree:split_bbox_if_flipped({160,-40,-120,60}, {-170,-80,170,80}),
            [{-170,-40,-120,60}, {160,-40,170,60}],
            "Bbox (flipped in x direction) with different bounds"),
    etap:is(vtree:split_bbox_if_flipped({-30,70,-20,-60}, {-170,-80,170,80}),
            [{-30,-80,-20,-60}, {-30,70,-20,80}],
            "Bbox (flipped in x direction) with different bounds"),
    etap:is(vtree:split_bbox_if_flipped({150,70,30,-60}, {-170,-80,170,80}),
            [{-170,-80,30,-60}, {-170,70,30,80},
             {150,-80,170,-60}, {150,70,170,80}],
            "Bbox (flipped in x and y direction) with different bounds"),

    etap:is(vtree:split_bbox_if_flipped(
             {-180, 28, 180, -28}, {-20, -20, 20, 20}), [],
           "Bbox that would be flipped in y direction," ++
           "but is out of bounds"),
    etap:is(vtree:split_bbox_if_flipped(
             {28, -90, -28, 90}, {-20, -20, 20, 20}), [],
           "Bbox that would be flipped in x direction," ++
           "but is out of bounds"),
    etap:is(vtree:split_bbox_if_flipped(
             {28, 28, -28, -28}, {-20, -20, 20, 20}), [],
           "Bbox that would be flipped in x and y direction," ++
           "but is out of bounds"),

    etap:is(vtree:split_bbox_if_flipped(
             {-180, 15, 180, -28}, {-20, -20, 20, 20}), [{-180,15,180,20}],
           "Bbox that is flipped in y direction," ++
           "one side is out of bounds"),
    etap:is(vtree:split_bbox_if_flipped(
             {28, -90, -18, 90}, {-20, -20, 20, 20}), [{-20,-90,-18,90}],
           "Bbox that would be flipped in x direction," ++
           "one side is out of bounds"),
    etap:is(vtree:split_bbox_if_flipped(
             {18, 28, -28, -15}, {-20, -20, 20, 20}), [{18,-20,20,-15}],
           "Bbox that would be flipped in x and y direction," ++
           "one side of each direction is out of bounds"),
    ok.


test_area() ->
    %etap:plan(5),
    Mbr1 = {10,5,13,15},
    Mbr2 = {-18,-3,-10,-1},
    Mbr3 = {-21,2,-10,14},
    Mbr4 = {5,-32,19,-25},
    Mbr5 = {-5,-16,4,19},
    etap:is(vtree:area(Mbr1), 30, "Area of MBR in the NE"),
    etap:is(vtree:area(Mbr2), 16, "Area of MBR in the SW"),
    etap:is(vtree:area(Mbr3), 132, "Area of MBR in the NW"),
    etap:is(vtree:area(Mbr4), 98, "Area of MBR in the SE"),
    etap:is(vtree:area(Mbr5), 315, "Area of MBR covering all quadrants"),
    ok.

test_merge_mbr() ->
    %etap:plan(7),
    Mbr1 = {10,5,13,15},
    Mbr2 = {-18,-3,-10,-1},
    Mbr3 = {-21,2,-10,14},
    Mbr4 = {5,-32,19,-25},
    etap:is(vtree:merge_mbr(Mbr1, Mbr2), {-18, -3, 13, 15},
            "Merge MBR of MBRs in NE and SW"),
    etap:is(vtree:merge_mbr(Mbr1, Mbr3), {-21, 2, 13, 15},
            "Merge MBR of MBRs in NE and NW"),
    etap:is(vtree:merge_mbr(Mbr1, Mbr4), {5, -32, 19, 15},
            "Merge MBR of MBRs in NE and SE"),
    etap:is(vtree:merge_mbr(Mbr2, Mbr3), {-21, -3, -10, 14},
            "Merge MBR of MBRs in SW and NW"),
    etap:is(vtree:merge_mbr(Mbr2, Mbr4), {-18, -32, 19, -1},
            "Merge MBR of MBRs in SW and SE"),
    etap:is(vtree:merge_mbr(Mbr3, Mbr4), {-21, -32, 19, 14},
            "Merge MBR of MBRs in NW and SE"),
    etap:is(vtree:merge_mbr(Mbr1, Mbr1), Mbr1,
            "Merge MBR of equal MBRs"),
    ok.

test_find_area_min_nth() ->
    %etap:plan(5),
    etap:is(vtree:find_area_min_nth([{5, {23,64,24,79}}]), 1,
            "Find position of minimum area in a list with one element"),
    etap:is(vtree:find_area_min_nth([{538, {2,64,4,79}}, {29, {2,64,4,79}}]), 2,
            "Find position of minimum area in a list with two elements (1>2)"),
    etap:is(vtree:find_area_min_nth([{54, {2,64,4,79}}, {538, {2,64,4,79}}]), 1,
            "Find position of minimum area in a list with two elements (1<2)"),
    etap:is(vtree:find_area_min_nth([{54, {2,64,4,79}}, {54, {2,64,4,79}}]), 1,
            "Find position of minimum area in a list with two equal elements"),
    etap:is(vtree:find_area_min_nth(
              [{329, {2,64,4,79}}, {930, {2,64,4,79}}, {203, {2,64,4,79}},
               {72, {2,64,4,79}}, {402, {2,64,4,79}}, {2904, {2,64,4,79}},
               {283, {2,64,4,79}}]), 4,
            "Find position of minimum area in a list"),
    ok.

test_partition_node() ->
    %etap:plan(3),
    Node1 = {{10,5,13,15}, #node{type=leaf}, <<"Node1">>},
    Node2 = {{-18,-3,-10,-1}, #node{type=leaf}, <<"Node2">>},
    Node3 = {{-21,2,-10,14}, #node{type=leaf}, <<"Node3">>},
    Node4 = {{5,-32,19,-25}, #node{type=leaf}, <<"Node4">>},
    Node5 = {{-5,-16,4,19}, #node{type=leaf}, <<"Node5">>},
    Children3 = [Node1, Node2, Node4],
    Children4 = Children3 ++ [Node3],
    Children5 = Children4 ++ [Node5],
    Mbr1_2_4 = {-18,-25,19,15},
    Mbr1_2_3_4_5 = {-21,-25,19,19},
    etap:is(vtree:partition_node({Mbr1_2_4, #node{type=leaf}, Children3}),
            {[Node2], [Node4], [Node1, Node4], [Node1, Node2]},
            "Partition 3 nodes"),
    etap:is(vtree:partition_node({Mbr1_2_3_4_5, #node{type=leaf}, Children4}),
            {[Node2, Node3], [Node4],
             [Node1, Node4], [Node1, Node2, Node3]},
            "Partition 4 nodes"),
    etap:is(vtree:partition_node({Mbr1_2_3_4_5, #node{type=leaf}, Children5}),
            {[Node2, Node3], [Node4],
             [Node1, Node4, Node5], [Node1, Node2, Node3, Node5]},
            "Partition 5 nodes"),
    ok.


test_calc_mbr() ->
    %etap:plan(10),
    Mbr1 = {10,5,13,15},
    Mbr2 = {-18,-3,-10,-1},
    Mbr3 = {-21,2,-10,14},
    Mbr4 = {5,-32,19,-25},
    etap:is(vtree:calc_mbr([]), error,
            "Calculate MBR of an empty list"),
    etap:is(vtree:calc_mbr([Mbr1]), {10, 5, 13, 15},
            "Calculate MBR of a single MBR"),
    etap:is(vtree:calc_mbr([Mbr1, Mbr2]), {-18, -3, 13, 15},
            "Calculate MBR of MBRs in NE and SW"),
    etap:is(vtree:calc_mbr([Mbr1, Mbr3]), {-21, 2, 13, 15},
            "Calculate MBR of MBRs in NE and NW"),
    etap:is(vtree:calc_mbr([Mbr1, Mbr4]), {5, -32, 19, 15},
            "Calculate MBR of MBRs in NE and SE"),
    etap:is(vtree:calc_mbr([Mbr2, Mbr3]), {-21, -3, -10, 14},
            "Calculate MBR of MBRs in SW and NW"),
    etap:is(vtree:calc_mbr([Mbr2, Mbr4]), {-18, -32, 19, -1},
            "Calculate MBR of MBRs in SW and SE"),
    etap:is(vtree:calc_mbr([Mbr3, Mbr4]), {-21, -32, 19, 14},
            "Calculate MBR of MBRs in NW and SE"),
    etap:is(vtree:calc_mbr([Mbr1, Mbr2, Mbr3]), {-21, -3, 13, 15},
            "Calculate MBR of MBRs in NE, SW, NW"),
    etap:is(vtree:calc_mbr([Mbr1, Mbr2, Mbr4]), {-18, -32, 19, 15},
            "Calculate MBR of MBRs in NE, SW, SE"),
    ok.

test_calc_nodes_mbr() ->
    %etap:plan(9),
    Node1 = {{10,5,13,15}, #node{type=leaf}, <<"Node1">>},
    Node2 = {{-18,-3,-10,-1}, #node{type=leaf}, <<"Node2">>},
    Node3 = {{-21,2,-10,14}, #node{type=leaf}, <<"Node3">>},
    Node4 = {{5,-32,19,-25}, #node{type=leaf}, <<"Node4">>},
    etap:is(vtree:calc_nodes_mbr([Node1]), {10, 5, 13, 15},
            "Calculate MBR of a single nodes"),
    etap:is(vtree:calc_nodes_mbr([Node1, Node2]), {-18, -3, 13, 15},
            "Calculate MBR of nodes in NE and SW"),
    etap:is(vtree:calc_nodes_mbr([Node1, Node3]), {-21, 2, 13, 15},
            "Calculate MBR of nodes in NE and NW"),
    etap:is(vtree:calc_nodes_mbr([Node1, Node4]), {5, -32, 19, 15},
            "Calculate MBR of nodes in NE and SE"),
    etap:is(vtree:calc_nodes_mbr([Node2, Node3]), {-21, -3, -10, 14},
            "Calculate MBR of nodes in SW and NW"),
    etap:is(vtree:calc_nodes_mbr([Node2, Node4]), {-18, -32, 19, -1},
            "Calculate MBR of nodes in SW and SE"),
    etap:is(vtree:calc_nodes_mbr([Node3, Node4]), {-21, -32, 19, 14},
            "Calculate MBR of nodes in NW and SE"),
    etap:is(vtree:calc_nodes_mbr([Node1, Node2, Node3]), {-21, -3, 13, 15},
            "Calculate MBR of nodes in NE, SW, NW"),
    etap:is(vtree:calc_nodes_mbr([Node1, Node2, Node4]), {-18, -32, 19, 15},
            "Calculate MBR of nodes in NE, SW, SE"),
    ok.

test_best_split() ->
    %etap:plan(4),
    Node1 = {{10,5,13,15}, #node{type=leaf}, <<"Node1">>},
    Node2 = {{-18,-3,-10,-1}, #node{type=leaf}, <<"Node2">>},
    Node3 = {{-21,2,-10,14}, #node{type=leaf}, <<"Node3">>},
    Node4 = {{5,-32,19,-25}, #node{type=leaf}, <<"Node4">>},
    Node5 = {{-5,-16,4,19}, #node{type=leaf}, <<"Node5">>},
    Node6 = {{-15,10,-12,17}, #node{type=leaf}, <<"Node6">>},
    {Mbr2, _, _} = Node2,
    {Mbr4, _, _} = Node4,
    Mbr1_2 = {-18,-3,13,15},
    Mbr1_4 = {5,-32,19,15},
    Mbr1_4_5 = {-5,-32,19,19},
    Mbr2_3 = {-21,-3,-10,14},
    Mbr4_6 = {-15,-32,19,17},
    Partition3 = {[Node2], [Node4], [Node1, Node4], [Node1, Node2]},
    Partition4 = {[Node2, Node3], [Node4],
                  [Node1, Node4], [Node1, Node2, Node3]},
    Partition5 = {[Node2, Node3], [Node4],
                  [Node1, Node4, Node5], [Node1, Node2, Node3, Node5]},
    Partition4b = {[Node2], [Node4, Node6],
                   [Node1, Node4, Node6], [Node1, Node2]},
    etap:is(vtree:best_split(Partition3), {tie, {Mbr2, Mbr4, Mbr1_4, Mbr1_2}},
            "Best split: tie (3 nodes)"),
    etap:is(vtree:best_split(Partition4), [{Mbr2_3, [Node2, Node3]},
                                           {Mbr1_4, [Node1, Node4]}],
            "Best split: horizontal (W/E) nodes win (4 nodes)"),
    etap:is(vtree:best_split(Partition5), [{Mbr2_3, [Node2, Node3]},
                                           {Mbr1_4_5, [Node1, Node4, Node5]}],
            "Best split: horizontal (W/E) nodes win (5 nodes)"),
    etap:is(vtree:best_split(Partition4b), [{Mbr4_6, [Node4, Node6]},
                                            {Mbr1_2, [Node1, Node2]}],
            "Best split: vertical (S/N) nodes win (4 nodes)"),
    ok.

test_minimal_overlap() ->
    % XXX vmx: test fir S/N split is missing
    %etap:plan(2),
    Node1 = {{10,5,13,15}, <<"Node1">>},
    Node2 = {{-18,-3,-10,-1}, <<"Node2">>},
    Node3 = {{-21,2,-10,14}, <<"Node3">>},
    Node4 = {{5,-32,19,-25}, <<"Node4">>},
    Node5 = {{-11,-9,12,10}, <<"Node5">>},
    {Mbr2, _} = Node2,
    {Mbr4, _} = Node4,
    Mbr1_2 = {-18,-3,13,15},
    Mbr1_3 = {-21,2,13,15},
    Mbr1_4 = {5,-32,19,15},
    Mbr1_4_5 = {-11,-32,19,15},
    Mbr2_3 = {-21,-3,-10,14},
    Mbr2_4_5 = {-18,-32,19,10},
    Partition3 = {[Node2], [Node4], [Node1, Node4], [Node1, Node2]},
    Partition5 = {[Node2, Node3], [Node2, Node4, Node5],
                  [Node1, Node4, Node5], [Node1, Node3]},
    etap:is(vtree:minimal_overlap(
                Partition5, {Mbr2_3, Mbr2_4_5, Mbr1_4_5, Mbr1_3}),
            [{Mbr2_3, [Node2, Node3]}, {Mbr1_4_5, [Node1, Node4, Node5]}],
            "Minimal Overlap: horizontal (W/E) nodes win (5 Nodes)"),
    etap:is(vtree:minimal_overlap(Partition3, {Mbr2, Mbr4, Mbr1_4, Mbr1_2}),
            tie, "Minimal Overlap: tie"),
    ok.


test_minimal_coverage() ->
    % XXX vmx: test for equal coverage is missing
    %etap:plan(2),
    Node1 = {{10,5,13,15}, <<"Node1">>},
    Node2 = {{-18,-3,-10,-1}, <<"Node2">>},
    Node3 = {{-21,2,-10,14}, <<"Node3">>},
    Node4 = {{5,-32,19,-25}, <<"Node4">>},
    Node5 = {{-11,-9,12,10}, <<"Node5">>},
    Node6 = {{-11,-9,12,24}, <<"Node6">>},
    {Mbr4, _} = Node4,
    Mbr1_3 = {-21,2,13,15},
    Mbr1_2_3_6 = {-21,-9,13,24},
    Mbr1_4_5 = {-11,-32,19,15},
    Mbr1_4_6 = {-11,-32,19,24},
    Mbr2_3 = {-21,-3,-10,14},
    Mbr2_4_5 = {-18,-32,19,10},
    Partition5 = {[Node2, Node3], [Node2, Node4, Node5],
                  [Node1, Node4, Node5], [Node1, Node3]},
    Partition6 = {[Node2, Node3], [Node4],
                  [Node1, Node4, Node6], [Node1, Node2, Node3, Node6]},
    etap:is(vtree:minimal_coverage(
                Partition5, {Mbr2_3, Mbr2_4_5, Mbr1_4_5, Mbr1_3}),
            [{Mbr2_3, [Node2, Node3]}, {Mbr1_4_5, [Node1, Node4, Node5]}],
            "Minimal Overlap: horizontal (W/E) nodes win)"),
    etap:is(vtree:minimal_coverage(
                Partition6, {Mbr2_3, Mbr4, Mbr1_4_6, Mbr1_2_3_6}),
            [{Mbr4, [Node4]}, {Mbr1_2_3_6, [Node1, Node2, Node3, Node6]}],
            "Minimal Overlap: vertical (S/N) nodes win"),
    ok.


test_calc_overlap() ->
    %etap:plan(7),
    Mbr1 = {10,5,13,15},
    Mbr2 = {-18,-3,-10,-1},
    Mbr3 = {-21,2,-10,14},
    Mbr4 = {5,-32,19,-25},
    Mbr5 = {-11,-9,12,10},
    Mbr6 = {-5,-6,4,9},
    Mbr7 = {4,-11,20,-3},
    etap:is(vtree:calc_overlap(Mbr1, Mbr5), {10, 5, 12, 10},
            "Calculate overlap of MBRs in NE and center"),
    etap:is(vtree:calc_overlap(Mbr2, Mbr5), {-11, -3, -10, -1},
            "Calculate overlap of MBRs in SW and center"),
    etap:is(vtree:calc_overlap(Mbr3, Mbr5), {-11, 2, -10, 10},
            "Calculate overlap of MBRs in NW and center"),
    etap:is(vtree:calc_overlap(Mbr7, Mbr5), {4, -9, 12, -3},
            "Calculate overlap of MBRs in SE and center"),
    etap:is(vtree:calc_overlap(Mbr6, Mbr5), {-5, -6, 4, 9},
            "Calculate overlap of one MBRs enclosing the other (1)"),
    etap:is(vtree:calc_overlap(Mbr5, Mbr6), {-5, -6, 4, 9},
            "Calculate overlap of one MBRs enclosing the other (2)"),
    etap:is(vtree:calc_overlap(Mbr4, Mbr5), {0, 0, 0, 0},
            "Calculate overlap of MBRs with no overlap"),
    ok.

test_delete() ->
    %etap:plan(9),

    {ok, Fd} = case couch_file:open(?FILENAME, [create, overwrite]) of
    {ok, Fd2} ->
        {ok, Fd2};
    {error, Reason} ->
        io:format("ERROR (~s): Couldn't open file (~s) for tree storage~n",
                  [Reason, ?FILENAME])
    end,

    Node1 = {{10,5,13,15}, #node{type=leaf},
        {linestring, [[10,5],[13,15]]}, <<"Node1">>},
    Node2 = {{-18,-3,-10,-1}, #node{type=leaf},
        {linestring, [[-18,-3],[-10,-1]]}, <<"Node2">>},
    Node3 = {{-21,2,-10,14}, #node{type=leaf},
        {linestring, [[-21,2],[-10,14]]}, <<"Node3">>},
    Node4 = {{5,-32,19,-25}, #node{type=leaf},
        {linestring, [[5,-32],[19,-25]]}, <<"Node4">>},
    Node5 = {{-5,-16,4,19}, #node{type=leaf},
        {linestring, [[-5,-16],[4,19]]}, <<"Node5">>},
    % od = on disk
    Node1od = {{10,5,13,15}, #node{type=leaf},
        {<<"Node1">>, {{linestring, [[10,5],[13,15]]}, <<"Node1">>}}},
    Node2od = {{-18,-3,-10,-1}, #node{type=leaf},
        {<<"Node2">>, {{linestring, [[-18,-3],[-10,-1]]}, <<"Node2">>}}},
    Node3od = {{-21,2,-10,14}, #node{type=leaf},
        {<<"Node3">>, {{linestring, [[-21,2],[-10,14]]}, <<"Node3">>}}},
    Node4od = {{5,-32,19,-25}, #node{type=leaf},
        {<<"Node4">>, {{linestring, [[5,-32],[19,-25]]}, <<"Node4">>}}},
    Node5od = {{-5,-16,4,19}, #node{type=leaf},
        {<<"Node5">>, {{linestring, [[-5,-16],[4,19]]}, <<"Node5">>}}},
    Mbr1 = {10,5,13,15},
    Mbr1_2 = {-18,-3,13,15},
    Mbr1_2_3 = {-21,-3,13,15},
    Mbr1_2_3_4 = {-21,-32,19,15},
    Mbr1_2_3_4_5 = {-21,-32,19,19},
    Mbr1_3_4 = {-21,-32,19,15},
    Mbr1_3_4_5 = {-21,-32,19,19},
    Mbr1_4_5 = {-5,-32,19,19},
    Mbr1_4 = {5,-32,19,15},
    {Mbr2, _, _, _} = Node2,
    Mbr2_3 = {-21,-3,-10,14},
    Mbr2_3_4_5 = {-21,-32,19,19},
    {Mbr3, _, _, _} = Node3,
    Mbr4_5 = {-5,-32,19,19},

    {ok, Mbr1, 0, 1} = vtree:insert(Fd, nil, <<"Node1">>, Node1),
    {ok, Mbr1_2, Pos2, 1} = vtree:insert(Fd, 0, <<"Node2">>, Node2),
    {ok, Mbr1_2_3, Pos3, 1} = vtree:insert(Fd, Pos2, <<"Node3">>, Node3),
    {ok, Mbr1_2_3_4, Pos4, 1} = vtree:insert(Fd, Pos3, <<"Node4">>, Node4),
    {ok, Mbr1_2_3_4_5, Pos5, 2} = vtree:insert(Fd, Pos4, <<"Node5">>, Node5),
    Tree3 = {Mbr3, #node{type=leaf}, [Node3od]},
    Tree2_3 = {Mbr2_3, #node{type=leaf}, [Node2od, Node3od]},
    Tree4_5 = {Mbr4_5, #node{type=leaf}, [Node4od, Node5od]},
    Tree1_4_5 = {Mbr1_4_5, #node{type=leaf}, [Node1od, Node4od, Node5od]},
    Tree2_3_4_5 = {Mbr2_3_4_5, #node{type=inner}, [Tree2_3, Tree4_5]},
    Tree1_3_4_5 = {Mbr1_3_4_5, #node{type=inner}, [Tree1_4_5, Tree3]},
    Tree1_4 = {Mbr1_4, #node{type=leaf}, [Node1od, Node4od]},
    Tree1_3_4 = {Mbr1_3_4, #node{type=inner}, [Tree3, Tree1_4]},

    {Node1Mbr, _,  _, Node1Id} = Node1,
    {Node2Mbr, _,  _, Node2Id} = Node2,
    {Node3Mbr, _,  _, Node3Id} = Node3,
    {Node4Mbr, _,  _, Node4Id} = Node4,
    {Node5Mbr, _,  _, Node5Id} = Node5,
    etap:is(vtree:delete(Fd, <<"bliblablubfoobar">>, Node1Mbr, Pos2),
            not_found,
            "Delete a node which ID's doesn't exist (tree height=1)"),

    {ok, Pos2_1} = vtree:delete(Fd, Node1Id, Node1Mbr, Pos2),
    etap:is(get_node(Fd, Pos2_1),
            {ok, {Mbr2, #node{type=leaf}, [Node2od]}},
            "Delete a node (tree height=1) (a)"),

    {ok, Pos2_2} = vtree:delete(Fd, Node2Id, Node2Mbr, Pos2),
    etap:is(get_node(Fd, Pos2_2),
            {ok, {Mbr1, #node{type=leaf}, [Node1od]}},
            "Delete a node (tree height=1) (b)"),

    {ok, Pos5_1} = vtree:delete(Fd, Node1Id, Node1Mbr, Pos5),
    {ok, {Pos5_1Mbr, Pos5_1Meta, [Pos5_1C1, Pos5_1C2]}} = get_node(
                                                            Fd, Pos5_1),
    {ok, Pos5_1Child1} = get_node(Fd, Pos5_1C1),
    {ok, Pos5_1Child2} = get_node(Fd, Pos5_1C2),
    etap:is({Pos5_1Mbr, Pos5_1Meta, [Pos5_1Child1, Pos5_1Child2]}, Tree2_3_4_5,
            "Delete a node (tree height=2) (a)"),

    {ok, Pos5_2} = vtree:delete(Fd, Node2Id, Node2Mbr, Pos5),
    {ok, {Pos5_2Mbr, Pos5_2Meta, [Pos5_2C1, Pos5_2C2]}} = get_node(
                                                            Fd, Pos5_2),
    {ok, Pos5_2Child1} = get_node(Fd, Pos5_2C1),
    {ok, Pos5_2Child2} = get_node(Fd, Pos5_2C2),
    etap:is({Pos5_2Mbr, Pos5_2Meta, [Pos5_2Child1, Pos5_2Child2]}, Tree1_3_4_5,
            "Delete a node (tree height=2) (b)"),

    {ok, Pos5_3} = vtree:delete(Fd, Node3Id, Node3Mbr, Pos5_2),
    {ok, {Pos5_3Mbr, Pos5_3Meta, [Pos5_3C]}} = get_node(Fd, Pos5_3),
    {ok, Pos5_3Child} = get_node(Fd, Pos5_3C),
    etap:is({Pos5_3Mbr, Pos5_3Meta, [Pos5_3Child]},
            {Mbr1_4_5, #node{type=inner}, [Tree1_4_5]},
            "Delete a node which is the only child (tree height=2) (b)"),

    {ok, Pos5_4} = vtree:delete(Fd, Node5Id, Node5Mbr, Pos5_2),
    {ok, {Pos5_4Mbr, Pos5_4Meta, [Pos5_4C1, Pos5_4C2]}} = get_node(
                                                            Fd, Pos5_4),
    {ok, Pos5_4Child1} = get_node(Fd, Pos5_4C1),
    {ok, Pos5_4Child2} = get_node(Fd, Pos5_4C2),
    etap:is({Pos5_4Mbr, Pos5_4Meta, [Pos5_4Child1, Pos5_4Child2]}, Tree1_3_4,
            "Delete a node (tree height=2) (b)"),

    % previous tests test the same code path already
    {ok, Pos5_5} = vtree:delete(Fd, Node4Id, Node4Mbr, Pos5_4),
    {ok, Pos5_6} = vtree:delete(Fd, Node3Id, Node3Mbr, Pos5_5),

    etap:is(vtree:delete(Fd, Node1Id, Node1Mbr, Pos5_6), {empty, nil},
            "After deletion of node, the tree is empty (tree height=2)"),

    etap:is(vtree:delete(Fd, Node5Id, Node5Mbr, Pos5_6), not_found,
            "Node can't be found (tree height=2)"),


    % Test what happens if all nodes have the same MBR
    Mbr = {-10,-10,20,20},
    Node11 = {Mbr, #node{type=leaf}, {linestring, [[-10,-10],[20,20]]},
        <<"Node11">>},
    Node12 = {Mbr, #node{type=leaf}, {linestring, [[-10,-10],[20,20]]},
        <<"Node12">>},
    Node13 = {Mbr, #node{type=leaf}, {linestring, [[-10,-10],[20,20]]},
        <<"Node13">>},
    Node14 = {Mbr, #node{type=leaf}, {linestring, [[-10,-10],[20,20]]},
        <<"Node14">>},
    Node15 = {Mbr, #node{type=leaf}, {linestring, [[-10,-10],[20,20]]},
        <<"Node15">>},
    % od = on disk
    Node11od = {Mbr, #node{type=leaf}, {<<"Node11">>,
        {{linestring, [[-10,-10],[20,20]]}, <<"Node11">>}}},
    Node12od = {Mbr, #node{type=leaf}, {<<"Node12">>,
        {{linestring, [[-10,-10],[20,20]]}, <<"Node12">>}}},
    Node13od = {Mbr, #node{type=leaf}, {<<"Node13">>,
        {{linestring, [[-10,-10],[20,20]]}, <<"Node13">>}}},
    Node14od = {Mbr, #node{type=leaf}, {<<"Node14">>,
        {{linestring, [[-10,-10],[20,20]]}, <<"Node14">>}}},
    {ok, Mbr, Pos11, 1} = vtree:insert(Fd, nil, <<"Node11">>, Node11),
    {ok, Mbr, Pos12, 1} = vtree:insert(Fd, Pos11, <<"Node12">>, Node12),
    {ok, Mbr, Pos13, 1} = vtree:insert(Fd, Pos12, <<"Node13">>, Node13),
    {ok, Mbr, Pos14, 1} = vtree:insert(Fd, Pos13, <<"Node14">>, Node14),
    {ok, Mbr, Pos15, 2} = vtree:insert(Fd, Pos14, <<"Node15">>, Node15),

    {ok, Pos16} = vtree:delete(Fd, <<"Node15">>, Mbr, Pos15),
    {ok, {Mbr, _, [Pos16C1, Pos16C2]}} = get_node(Fd, Pos16),
    {ok, {Mbr, _, [Node11od, Node12od]}} = get_node(Fd, Pos16C1),
    {ok, {Mbr, _, [Node13od, Node14od]}} = get_node(Fd, Pos16C2),


    {ok, Pos17} = vtree:delete(Fd, <<"Node12">>, Mbr, Pos16),
    {ok, {Mbr, _, [Pos17C1, Pos17C2]}} = get_node(Fd, Pos17),
    {ok, {Mbr, _, [Node13od, Node14od]}} = get_node(Fd, Pos17C1),
    {ok, {Mbr, _, [Node11od]}} = get_node(Fd, Pos17C2),

    {ok, Pos18} = vtree:delete(Fd, <<"Node11">>, Mbr, Pos17),
    {ok, {Mbr, _, [Pos18C1]}} = get_node(Fd, Pos18),
    {ok, {Mbr, _, [Node13od, Node14od]}} = get_node(Fd, Pos18C1),

    {ok, Pos19} = vtree:delete(Fd, <<"Node13">>, Mbr, Pos18),
    {ok, {Mbr, _, [Pos19C1]}} = get_node(Fd, Pos19),
    {ok, {Mbr, _, [Node14od]}} = get_node(Fd, Pos19C1),

    {empty, nil} = vtree:delete(Fd, <<"Node14">>, Mbr, Pos19),
    ok.

test_delete_same_id() ->
    % Test what happens if multiple emits in one function happened
    %etap:plan(9),

    {ok, Fd} = case couch_file:open(?FILENAME, [create, overwrite]) of
    {ok, Fd2} ->
        {ok, Fd2};
    {error, Reason} ->
        io:format("ERROR (~s): Couldn't open file (~s) for tree storage~n",
                  [Reason, ?FILENAME])
    end,

    Node21 = {{10,5,13,15}, #node{type=leaf},
        {linestring, [[10,5],[13,15]]}, <<"Node">>},
    Node22 = {{-18,-3,-10,-1}, #node{type=leaf},
        {linestring, [[-18,-3],[-10,-1]]}, <<"Node">>},
    Node23 = {{-21,2,-10,14}, #node{type=leaf},
        {linestring, [[-21,2],[-10,14]]}, <<"Node">>},
    Node24 = {{5,-32,19,-25}, #node{type=leaf},
        {linestring, [[5,-32],[19,-25]]}, <<"Node">>},
    Node25 = {{-5,-16,4,19}, #node{type=leaf},
        {linestring, [[-5,-16],[4,19]]}, <<"Node">>},
    % od = on disk
    Node21od = {{10,5,13,15}, #node{type=leaf},
        {<<"Node">>, {{linestring, [[10,5],[13,15]]}, <<"Node">>}}},
    Node22od = {{-18,-3,-10,-1}, #node{type=leaf},
        {<<"Node">>, {{linestring, [[-18,-3],[-10,-1]]}, <<"Node">>}}},
    Node23od = {{-21,2,-10,14}, #node{type=leaf},
        {<<"Node">>, {{linestring, [[-21,2],[-10,14]]}, <<"Node">>}}},
    Node24od = {{5,-32,19,-25}, #node{type=leaf},
        {<<"Node">>, {{linestring, [[5,-32],[19,-25]]}, <<"Node">>}}},
    Node25od = {{-5,-16,4,19}, #node{type=leaf},
        {<<"Node">>, {{linestring, [[-5,-16],[4,19]]}, <<"Node">>}}},

    Mbr1 = {10,5,13,15},
    Mbr1_2 = {-18,-3,13,15},
    Mbr1_2_3 = {-21,-3,13,15},
    Mbr1_2_3_4 = {-21,-32,19,15},
    Mbr1_2_3_4_5 = {-21,-32,19,19},
    Mbr1_3_4 = {-21,-32,19,15},
    Mbr1_3_4_5 = {-21,-32,19,19},
    Mbr1_3_4 = {-21,-32,19,15},
    Mbr1_3_4_5 = {-21,-32,19,19},
    Mbr1_4_5 = {-5,-32,19,19},
    Mbr1_4 = {5,-32,19,15},
    {Mbr2, _, _, _} = Node22,
    Mbr2_3 = {-21,-3,-10,14},
    Mbr2_3_4_5 = {-21,-32,19,19},
    {Mbr3, _, _, _} = Node23,
    Mbr4_5 = {-5,-32,19,19},

    {ok, Mbr1, 0, 1} = vtree:insert(Fd, nil, <<"Node">>, Node21),
    {ok, Mbr1_2, Pos22, 1} = vtree:insert(Fd, 0, <<"Node">>, Node22),
    {ok, Mbr1_2_3, Pos23, 1} = vtree:insert(Fd, Pos22, <<"Node">>, Node23),
    {ok, Mbr1_2_3_4, Pos24, 1} = vtree:insert(Fd, Pos23, <<"Node">>, Node24),
    {ok, Mbr1_2_3_4_5, Pos25, 2} = vtree:insert(Fd, Pos24, <<"Node">>, Node25),
    Tree23 = {Mbr3, #node{type=leaf}, [Node23od]},
    Tree22_3 = {Mbr2_3, #node{type=leaf}, [Node22od, Node23od]},
    Tree24_5 = {Mbr4_5, #node{type=leaf}, [Node24od, Node25od]},
    Tree21_4_5 = {Mbr1_4_5, #node{type=leaf}, [Node21od, Node24od, Node25od]},
    Tree22_3_4_5 = {Mbr2_3_4_5, #node{type=inner}, [Tree22_3, Tree24_5]},
    Tree21_3_4_5 = {Mbr1_3_4_5, #node{type=inner}, [Tree21_4_5, Tree23]},
    Tree21_4 = {Mbr1_4, #node{type=leaf}, [Node21od, Node24od]},
    Tree21_3_4 = {Mbr1_3_4, #node{type=inner}, [Tree23, Tree21_4]},

    {Node1Mbr, _, _, <<"Node">>} = Node21,
    {Node2Mbr, _, _, <<"Node">>} = Node22,
    {Node3Mbr, _, _, <<"Node">>} = Node23,
    {Node4Mbr, _, _, <<"Node">>} = Node24,
    {Node5Mbr, _, _, <<"Node">>} = Node25,

    etap:is(vtree:delete(Fd, <<"bliblablubfoobar">>, Node1Mbr, Pos22),
            not_found,
            "Delete a node which ID's doesn't exist (tree height=1) " ++
            "(same ID)"),
    {ok, Pos22_1} = vtree:delete(Fd, <<"Node">>, Node1Mbr, Pos22),
    etap:is(get_node(Fd, Pos22_1),
            {ok, {Mbr2, #node{type=leaf}, [Node22od]}},
            "Delete a node (tree height=1) (a) (same ID)"),

    {ok, Pos22_2} = vtree:delete(Fd, <<"Node">>, Node2Mbr, Pos22),
    etap:is(get_node(Fd, Pos22_2),
            {ok, {Mbr1, #node{type=leaf}, [Node21od]}},
            "Delete a node (tree height=1) (b) (same ID)"),

    {ok, Pos25_1} = vtree:delete(Fd, <<"Node">>, Node1Mbr, Pos25),
    {ok, {Pos25_1Mbr, Pos25_1Meta, [Pos25_1C1, Pos25_1C2]}} = get_node(
                                                            Fd, Pos25_1),
    {ok, Pos25_1Child1} = get_node(Fd, Pos25_1C1),
    {ok, Pos25_1Child2} = get_node(Fd, Pos25_1C2),
    etap:is({Pos25_1Mbr, Pos25_1Meta, [Pos25_1Child1, Pos25_1Child2]},
            Tree22_3_4_5, "Delete a node (tree height=2) (a) (same ID)"),

    {ok, Pos25_2} = vtree:delete(Fd, <<"Node">>, Node2Mbr, Pos25),
    {ok, {Pos25_2Mbr, Pos25_2Meta, [Pos25_2C1, Pos25_2C2]}} = get_node(
                                                            Fd, Pos25_2),
    {ok, Pos25_2Child1} = get_node(Fd, Pos25_2C1),
    {ok, Pos25_2Child2} = get_node(Fd, Pos25_2C2),
    etap:is({Pos25_2Mbr, Pos25_2Meta, [Pos25_2Child1, Pos25_2Child2]},
            Tree21_3_4_5, "Delete a node (tree height=2) (b) (same ID)"),

    {ok, Pos25_3} = vtree:delete(Fd, <<"Node">>, Node3Mbr, Pos25_2),
    {ok, {Pos25_3Mbr, Pos25_3Meta, [Pos25_3C]}} = get_node(Fd, Pos25_3),
    {ok, Pos25_3Child} = get_node(Fd, Pos25_3C),
    etap:is({Pos25_3Mbr, Pos25_3Meta, [Pos25_3Child]},
            {Mbr1_4_5, #node{type=inner}, [Tree21_4_5]},
            "Delete a node which is the only child (tree height=2) (b) " ++
            "(same ID)"),

    {ok, Pos25_4} = vtree:delete(Fd, <<"Node">>, Node5Mbr, Pos25_2),
    {ok, {Pos25_4Mbr, Pos25_4Meta, [Pos25_4C1, Pos25_4C2]}} = get_node(
                                                            Fd, Pos25_4),
    {ok, Pos25_4Child1} = get_node(Fd, Pos25_4C1),
    {ok, Pos25_4Child2} = get_node(Fd, Pos25_4C2),
    etap:is({Pos25_4Mbr, Pos25_4Meta, [Pos25_4Child1, Pos25_4Child2]},
            Tree21_3_4, "Delete a node (tree height=2) (b) (same ID)"),

    % previous tests test the same code path already
    {ok, Pos25_5} = vtree:delete(Fd, <<"Node">>, Node4Mbr, Pos25_4),
    {ok, Pos25_6} = vtree:delete(Fd, <<"Node">>, Node3Mbr, Pos25_5),

    etap:is(vtree:delete(Fd, <<"Node">>, Node1Mbr, Pos25_6), {empty, nil},
            "After deletion of node, the tree is empty (tree height=2) " ++
            "(same ID)"),

    etap:is(vtree:delete(Fd, <<"Node">>, Node5Mbr, Pos25_6), not_found,
            "Node can't be found (tree height=2) (same ID)"),
    ok.

test_split_node() ->
    %etap:plan(3),

    NodeToSplit1 =  {{52,45,969,960},{node,leaf}, [
        {{52,320,597,856},{node,leaf},{<<"Node7">>,<<"Node7">>}},
        {{270,179,584,331},{node,leaf},{<<"Node9">>,<<"Node9">>}},
        {{441,502,580,960},{node,leaf},{<<"Node13">>,<<"Node13">>}},
        {{462,520,543,911},{node,leaf},{<<"Node15">>,<<"Node15">>}},
        {{493,45,969,938},{node,leaf},{<<"Node17">>,<<"Node17">>}}]},

    {SplittedMbr1, Node1a, Node1b} = vtree:split_node(NodeToSplit1),
    etap:is(SplittedMbr1, vtree:calc_nodes_mbr(element(3, NodeToSplit1)),
        "Split node with one item too much (MBR is right)"),
    etap:is(length(element(3, Node1a)), 3,
        "Split node with one item too much (correct number of children (a))"),
    etap:is(length(element(3, Node1b)), 2,
        "Split node with one item too much (correct number of children (b))"),
    ok.

test_count_total() ->
    %etap:plan(3),

    {ok, {Fd1, {RootPos1, _}}} = gc_test_util:build_random_tree(
        "/tmp/randtree.bin", 20),
    Count1 = vtree:count_total(Fd1, RootPos1),
    etap:is(Count1, 20, "Total number of geometries is correct (a)"),

    {ok, {Fd2, {RootPos2, _}}} = gc_test_util:build_random_tree(
        "/tmp/randtree.bin", 338),
    Count2 = vtree:count_total(Fd2, RootPos2),
    etap:is(Count2, 338, "Total number of geometries is correct (b)"),

    Count3 = vtree:count_total(Fd2, nil),
    etap:is(Count3, 0,
            "Total number of geometries is correct (for empty tree)"),
    ok.


%% Helpers %%

get_node(Fd, Pos) ->
    couch_file:pread_term(Fd, Pos).
