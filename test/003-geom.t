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

% tolerance for `almost_equal` comparisons
-define(TOLERANCE, 0.000000000001).
% mean earth radius in km
-define(R, 6370.986).

main(_) ->
    code:add_pathz(filename:dirname(escript:script_name())),
    gc_test_util:init_code_path(),
    etap:plan(41),
    case (catch test()) of
        ok ->
            etap:end_tests();
        Other ->
            etap:diag(io_lib:format("Test died abnormally: ~p", [Other])),
            etap:bail(Other)
    end,
    ok.

test() ->
    test_distance(),
    test_sphere_distance(),
    test_sphere_distance_points(),
    test_to_geographic_point(),
    test_within(),
    test_remainder(),

    %etap:end_tests().
    ok.


test_distance() ->
    %etap:plan(10),

    etap:is(geom:distance({5,5}, {0,0,10,10}, nil),
            0,
            "Point is inside Mbr"),

    etap:is(geom:distance({0,5}, {0,0,10,10}, nil),
            0,
            "Point is on the border of the Mbr"),

    etap:is(geom:distance({5,11}, {0,0,10,10}, nil),
            1,
            "Point is over the Mbr"),

    etap:is(geom:distance({11,-1}, {0,0,10,10}, nil),
            geom:distance({11,-1}, {12,-12,22,-2}, nil),
            "Point has the same distance to two Mbrs"),

    etap:ok(geom:distance({-170,75}, {160,70,170,80}, {-180,-90,180,90}) <
            geom:distance({-170,75}, {160,70,170,80}, nil),
            "Bounds are used"),

    etap:is(geom:distance({1,3}, {5,2,6,4}, {0,0,6,4}) , 1,
            "Bounds: x + width"),
    etap:is(geom:distance({5.5,1}, {5,3,6,4}, {0,0,6,4}) , 1,
            "Bounds: y + height"),
    etap:is(geom:distance({1,0}, {5,3,6,4}, {0,0,6,4}) , 1,
            "Bounds: x + width and y + height"),
    etap:is(geom:distance({5.5,3}, {5,0,6,1}, {0,0,6,4}) , 1,
            "Bounds: y - height"),
    etap:is(geom:distance({0,3}, {5,0,6,1}, {0,0,6,4}) , 1,
            "Bounds: x + width and y - height"),
    ok.


test_sphere_distance() ->
    %etap:plan(8),

    almost_equal(
        geom:sphere_distance({5.5, 49.5}, {5, 49, 6, 50}),
        0.0,
        "point is inside"),
    almost_equal(
        geom:sphere_distance({6, 49}, {5, 49, 6, 50}),
        0.0,
        "point is one of the corner points"),
    almost_equal(
        geom:sphere_distance({4.5, 50.5}, {5, 49, 6, 50}) * ?R,
        geom:sphere_distance(
            geom:to_geographic_point({4.5, 50.5}), geom:to_geographic_point({5, 50})) * ?R,
        "point is closest to the upper-left corner"),
    almost_equal(geom:sphere_distance(
            {5.5, 50.25}, {5, 49, 6, 50}) * ?R, 
        27.679217608354296,
        "point is closest to a point on the upper edge"),
    almost_equal(geom:sphere_distance(
            {6.25, 49.5}, {5, 49, 6, 50}) * ?R, 
        18.05375922675459,
        "point is closest to a point on the right edge"),
    almost_equal(geom:sphere_distance(
            {5.5, 48.75}, {5, 49, 6, 50}) * ?R, 
        27.918785932650948,
        "point is closest to a point on the lower edge"),
    almost_equal(geom:sphere_distance(
            {4.75, 49.5}, {5, 49, 6, 50}) * ?R, 
        18.05375922675459,
        "point is closest to a point on the left edge"),
    almost_equal(geom:sphere_distance(
            {90, 1.1}, {-1, -1, 1, 1}) * ?R, 
        9894.229276213218,
        "point is near great circle through upper edge, but not on it"),

    ok.

test_sphere_distance_points() ->
    %etap:plan(5),

    almost_equal(geom:sphere_distance(
            geom:to_geographic_point({5, 5}), geom:to_geographic_point({5, 5})), 
        0.0,
        "distance is 0 when the from/to point is the same"),
    almost_equal(geom:sphere_distance(
            geom:to_geographic_point({6.37, 49.56}), geom:to_geographic_point({6.4, 49.2})) * ?R, 
        40.088953975074745,
        "small distance"),
    almost_equal(geom:sphere_distance(
            geom:to_geographic_point({6.37, 49.56}), geom:to_geographic_point({6.371, 49.562})) * ?R, 
        0.23379277325477713,
        "very small distance"),
    almost_equal(geom:sphere_distance(
            geom:to_geographic_point({-89, 45}), geom:to_geographic_point({90, -23})) * ?R, 
        17567.074360903247,
        "huge distance"),
    almost_equal(geom:sphere_distance(
            geom:to_geographic_point({160, 40}), geom:to_geographic_point({-130, 15})) * ?R, 
        7249.972591631049,
        "shortest path across the dateline"),

    ok.


test_to_geographic_point() ->
    %etap:plan(5),

    PI = math:pi(),
    PI_2 = PI / 2.0,

    almost_equal(geom:to_geographic_point({0, 0}), {0, 0}, ""),
    almost_equal(geom:to_geographic_point({-180, 90}), {PI, PI_2}, ""),
    almost_equal(geom:to_geographic_point({90, -45}), {PI_2, PI_2 / -2.0}, ""),
    almost_equal(geom:to_geographic_point({-90, 45}), {-1 * PI_2, PI_2 / 2.0}, ""),
    almost_equal(geom:to_geographic_point({-450, 405}), {-1 * PI_2, PI_2 / 2.0},
        "coordinates are nomalized"),

    ok.


test_within() ->
    %etap:plan(4),

    Bbox1 = {-20, -10, 30, 21},
    Bbox2 = {-20, -10, 0, 0},
    Mbr1_2 = {-18,-3,13,15},
    Node1 = {{10,5,13,15}, <<"Node1">>},
    {Node1Mbr, _} = Node1,
    etap:is(geom:within(Node1Mbr, Bbox1), true, "MBR is within the BBox"),
    etap:is(geom:within(Node1Mbr, Node1Mbr), true, "MBR is within itself"),
    etap:is(geom:within(Node1Mbr, Bbox2), false,
            "MBR is not at all within BBox"),
    etap:is(geom:within(Mbr1_2, Bbox2), false, "MBR intersects BBox"),

    ok.


test_remainder() ->
    %etap:plan(9),

    almost_equal(geom:remainder(3.5, 2.3), -1.1, "3.5 % 2.3 = -1.1"),
    almost_equal(geom:remainder(2.7, 2.7), 0.0, "2.7 % 2.7 = 0.0"),
    almost_equal(geom:remainder(4.9, 1.6), 0.1, "4.9 % 1.6 = 0.1"),
    almost_equal(geom:remainder(18.14, 2.17), 0.78, "18.14 % 2.17 = 0.78"),
    almost_equal(geom:remainder(5.0, 2.0), 1.0, "5.0 % 2.0 = 1.0"),
    almost_equal(geom:remainder(1.3, -1.4), -0.1, "1.3 % -1.4 = -0.1"),
    almost_equal(geom:remainder(7.0, 2.0), -1.0, "7.0 % 2.0 = -1.0"),
    almost_equal(geom:remainder(3.0, 2.0), -1.0, "3.0 % 2.0 = -1.0"),
    almost_equal(geom:remainder(3.0, -2.0), -1.0, "3.0 % -2.0 = -1.0"),

    ok.


almost_equal(Got = {X1, Y1}, Expected = {X2, Y2}, Desc) ->
    case ((abs(X1 - X2) =< ?TOLERANCE) and (abs(Y1 - Y2) =< ?TOLERANCE)) of
    true ->
        etap:ok(true, Desc);
    false ->
        etap:is(Got, Expected, Desc)
    end;

almost_equal(Got, Expected, Desc) ->
    case (abs(Got - Expected) =< ?TOLERANCE) of
    true ->
        etap:ok(true, Desc);
    false ->
        etap:is(Got, Expected, Desc)
    end.
