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

-module(geom).

% Most of the functions in this module for calculating the spherical distance
% were shamelessly "borrowed" from PostGIS (lwgeodetic.c, rev. 9128).

-export([distance/3, sphere_distance/2, within/2]).

-ifdef(makecheck).
-export([to_geographic_point/1, remainder/2]).
-endif.

-import(math, [cos/1, sin/1, pow/2, sqrt/1, atan2/2, asin/1]).

-define(PI, 3.141592653589793).
-define(PI_2, 1.5707963267948966).
-define(MAX_DISTANCE, 1000000000.0).
-define(TOLERANCE, 1.0e-12).


% Calculates the minimum cartesian distance (called MINDIST) between a point and
% a MBR, see "Nearest Neighbour Queries" by Roussopoulos et al.
% (http://www.cs.ucr.edu/~tsotras/cs236/F11/roussopoulosNN95.pdf) (def. 3, pg. 3).
distance(
    Point = {X, Y},
    Mbr,
    Bounds) ->

    case within(Point, Mbr) of
    true -> 0;
    false ->
        case Bounds == nil of
        true ->
            distance2(Point, Mbr);
        false ->
            % if bounds are given, also calculate the distance for 8 additional
            % points and return the minimum
            {XMinB, YMinB, XMaxB, YMaxB} = Bounds,
            Width = XMaxB - XMinB,
            Height = YMaxB - YMinB,

            lists:min(lists:map(
                fun({DeltaX, DeltaY}) ->
                    distance2({X + DeltaX, Y + DeltaY}, Mbr)
                end,
                [{0, 0}, {Width, 0}, {-Width, 0},
                 {0, Height}, {Width, Height}, {-Width, Height},
                 {0, -Height}, {Width, -Height}, {-Width, -Height}]))
        end
    end.


distance2({X, Y}, {XMin, YMin, XMax, YMax}) ->
    DistX = math:pow(abs(X - getR(X, XMin, XMax)), 2),
    DistY = math:pow(abs(Y - getR(Y, YMin, YMax)), 2),

    DistX + DistY.


getR(P, S, T) ->
    if
        P < S -> S;
        P > T -> T;
        true -> P
    end.


% Calculates the minimum spherical distance in radians between a point 
% and a minimum-bounding-rectangle.
%
% The coordinates are expected to be in WGS 84.
sphere_distance(
    Point, 
    Mbr = {XMin, YMin, XMax, YMax}) ->

    case within(Point, Mbr) of
    true -> 0;
    false ->
        % convert to 'geographic points'
        PointGP = to_geographic_point(Point),
        LowerLeftGP = to_geographic_point({XMin, YMin}),
        UpperLeftGP = to_geographic_point({XMin, YMax}),
        UpperRightGP = to_geographic_point({XMax, YMax}),
        LowerRightGP = to_geographic_point({XMax, YMin}),

        Distances = [
            % distances to the corner points
            sphere_distance(PointGP, LowerLeftGP),
            sphere_distance(PointGP, UpperLeftGP),
            sphere_distance(PointGP, UpperRightGP),
            sphere_distance(PointGP, LowerRightGP),

            % distances to the edges between the points
            sphere_distance(PointGP, {LowerLeftGP, UpperLeftGP}),
            sphere_distance(PointGP, {UpperLeftGP, UpperRightGP}),
            sphere_distance(PointGP, {UpperRightGP, LowerRightGP}),
            sphere_distance(PointGP, {LowerRightGP, LowerLeftGP})
        ],
        lists:min(Distances)
    end;

% Calculates the minimum distance in radians between a point and an edge.
%
% This function was adapted from `edge_distance_to_point(..)` (lwgeodetic.c).
% The main difference is that the distance to the edge vertices is not checked, only
% the distance to the point ON the edge that is closest to the given point.
%
% The coordinates are expected to be spherical coordinates (units of radians).
sphere_distance(
    PointGP, 
    {PointFromGP = {_, _}, PointToGP = {_, _}}) ->

    case geographic_point_equals(PointFromGP, PointToGP) of
    true -> ?MAX_DISTANCE;
    false ->
        % find a point on the great circle through `PointFromGP` and 
        % `PointToGP` that is closest to `PointGP`
        N1 = normalize(robust_cross_product(PointFromGP, PointToGP)),
        Point = geog2cart(PointGP),
        N2 = vector_scale(N1, dot_product(Point, N1)),
        NearestGP = cart2geog(normalize(vector_difference(Point, N2))),

        % now check that the point is really on the circle segment
        % between `PointFromGP` and `PointToGP`
        case edge_contains_point({PointFromGP, PointToGP}, NearestGP) of
            true -> sphere_distance(PointGP, NearestGP);
            false -> ?MAX_DISTANCE
        end
    end;

% Calculates the distance in radians between two points using the Havesine formula.
% Multiply by 6370.986 (mean earth radius) to get the distance in kilometers.
%
% The coordinates are expected to be spherical coordinates (units of radians).
sphere_distance({Lon1, Lat1}, {Lon2, Lat2}) ->

    D_Lon = Lon2 - Lon1,
    Cos_D_Lon = cos(D_Lon),
    Cos_Lat2 = cos(Lat2),
    Sin_Lat2 = sin(Lat2),
    Cos_Lat1 = cos(Lat1),
    Sin_Lat1 = sin(Lat1),

    A1 = pow(Cos_Lat2 * sin(D_Lon), 2),
    A2 = pow(Cos_Lat1 * Sin_Lat2 - Sin_Lat1 * Cos_Lat2 * Cos_D_Lon, 2),
    A = sqrt(A1 + A2),
    B = Sin_Lat1 * Sin_Lat2 + Cos_Lat1 * Cos_Lat2 * Cos_D_Lon,

    atan2(A, B).


% Converts the coordinates of a point to radian values,
% see geographic_point_init(lon, lat)
to_geographic_point({X, Y}) ->
    Lon = longitude_radians_normalize(deg2rad(X)),
    Lat = latitude_radians_normalize(deg2rad(Y)),

    {Lon, Lat}.


deg2rad(Degree) ->
    ?PI * Degree / 180.0.


% Convert a longitude to the range of -pi,pi
longitude_radians_normalize(Lon)
    when (Lon == -1.0 * ?PI) -> ?PI;
longitude_radians_normalize(Lon)
    when (Lon == -2.0 * ?PI) -> 0.0;
longitude_radians_normalize(Lon) ->

    Lon1 = case (Lon > 2.0 * ?PI) of
    true -> remainder(Lon, 2.0 * ?PI);
    false -> Lon
    end,

    Lon2 = case (Lon1 < -2.0 * ?PI) of
    true -> remainder(Lon1, -2.0 * ?PI);
    false -> Lon1
    end,

    Lon3 = case (Lon2 > ?PI) of
    true -> -2.0 * ?PI + Lon2;
    false -> Lon2
    end,

    Lon4 = case (Lon3 < -1.0 * ?PI) of
    true -> 2.0 * ?PI + Lon3;
    false -> Lon3
    end,

    Lon4.


% Convert a latitude to the range of -pi/2,pi/2
latitude_radians_normalize(Lat) ->

    Lat1 = case (Lat > 2.0 *  ?PI) of
    true -> remainder(Lat, 2.0 * ?PI);
    false -> Lat
    end,

    Lat2 = case (Lat1 < -2.0 *  ?PI) of
    true -> remainder(Lat1, -2.0 * ?PI);
    false -> Lat1
    end,

    Lat3 = case (Lat2 > ?PI) of
    true -> ?PI - Lat2;
    false -> Lat2
    end,

    Lat4 = case (Lat3 < -1.0 * ?PI) of
    true -> -1.0 * ?PI - Lat3;
    false -> Lat3
    end,

    Lat5 = case (Lat4 > ?PI_2) of
    true -> ?PI - Lat4;
    false -> Lat4
    end,

    Lat6 = case (Lat5 < -1.0 * ?PI_2) of
    true -> -1.0 * ?PI - Lat5;
    false -> Lat5
    end,

    Lat6.


geographic_point_equals({Lon1, Lat1}, {Lon2, Lat2}) ->
    (abs(Lon1 - Lon2) =< ?TOLERANCE) and (abs(Lat1 - Lat2) =< ?TOLERANCE).


% Computes the cross product of two vectors using their lat, lng representations.
% Good even for small distances between p and q.
robust_cross_product({Lon1, Lat1}, {Lon2, Lat2}) ->
    Lon_Qpp = (Lon2 + Lon1) / -2.0,
    Lon_Qmp = (Lon2 - Lon1) / 2.0,
    Sin_Lat1_Minus_Lat2 = sin(Lat1-Lat2),
    Sin_Lat1_Plus_Lat2 = sin(Lat1+Lat2),
    Sin_Lon_Qpp = sin(Lon_Qpp),
    Sin_Lon_Qmp = sin(Lon_Qmp),
    Cos_Lon_Qpp = cos(Lon_Qpp),
    Cos_Lon_Qmp = cos(Lon_Qmp),
    X = Sin_Lat1_Minus_Lat2 * Sin_Lon_Qpp * Cos_Lon_Qmp -
           Sin_Lat1_Plus_Lat2 * Cos_Lon_Qpp * Sin_Lon_Qmp,
    Y = Sin_Lat1_Minus_Lat2 * Cos_Lon_Qpp * Cos_Lon_Qmp +
           Sin_Lat1_Plus_Lat2 * Sin_Lon_Qpp * Sin_Lon_Qmp,
    Z = cos(Lat1) * cos(Lat2) * sin(Lon2-Lon1),
    {X, Y, Z}.


% Normalize to a unit vector.
normalize({X, Y, Z}) ->
    D = sqrt(X*X + Y*Y + Z*Z),
    case (abs(D) =< ?TOLERANCE) of
    true -> {0.0, 0.0, 0.0};
    false -> {X / D, Y / D, Z / D}
    end.


% Convert spherical coordinates to cartesion coordinates on unit sphere
geog2cart({Lon, Lat}) ->
    X = cos(Lat) * cos(Lon),
    Y = cos(Lat) * sin(Lon),
    Z = sin(Lat),
    {X, Y, Z}.


% Convert cartesion coordinates to spherical coordinates on unit sphere
cart2geog({X, Y, Z}) ->
    {atan2(Y, X), asin(Z)}.


% Calculate the dot product of two unit vectors
% (-1 == opposite, 0 == orthogonal, 1 == identical)
dot_product({X1, Y1, Z1}, {X2, Y2, Z2}) ->
    ((X1*X2) + (Y1*Y2) + (Z1*Z2)).


% Scale a vector out by a factor
vector_scale({X, Y, Z}, Scale) ->
    {X * Scale, Y * Scale, Z * Scale}.


% Calculate the sum of two vectors
vector_sum({X1, Y1, Z1}, {X2, Y2, Z2}) ->
    {X1 + X2, Y1 + Y2, Z1 + Z2}.


% Calculate the difference of two vectors
vector_difference({X1, Y1, Z1}, {X2, Y2, Z2}) ->
    {X1 - X2, Y1 - Y2, Z1 - Z2}.


% Returns true if the point `PointGP` is on the minor edge defined by the
% end points of `EdgeGP`.
edge_contains_point(EdgeGP, PointGP) ->
    (edge_point_in_cone(EdgeGP, PointGP) and edge_point_on_plane(EdgeGP, PointGP)).


% Returns true if the point PointGP is inside the cone defined by the
% two ends PointFromGP and PointToGP.
edge_point_in_cone({PointFromGP, PointToGP}, PointGP) ->
    Vs = {VsX, VsY, VsZ} = geog2cart(PointFromGP),
    Ve = {VeX, VeY, VeZ} = geog2cart(PointToGP),

    % Antipodal case, everything is inside.
    case ((VsX == (-1.0 * VeX)) and (VsY == (-1.0 * VeY)) and (VsZ == (-1.0 * VeZ))) of
    true -> true;
    false ->
        Vp = geog2cart(PointGP),
        % The normalized sum bisects the angle between start and end.
        Vcp = normalize(vector_sum(Vs, Ve)),

        % The projection of start onto the center defines the minimum similarity
        Vs_Dot_Vcp = dot_product(Vs, Vcp),
        % The projection of candidate p onto the center
        Vp_Dot_Vcp = dot_product(Vp, Vcp),

        % If PointGP is more similar than start then PointGP is inside the cone

        % We want to test that vp_dot_vcp is >= vs_dot_vcp but there are
        % numerical stability issues for values that are very very nearly
        % equal. Unfortunately there are also values of vp_dot_vcp that are legitimately
        % very close to but still less than vs_dot_vcp which we also need to catch.
        % The tolerance of 10-17 seems to do the trick on 32-bit and 64-bit architectures,
        % for the test cases here.
        % However, tuning the tolerance value feels like a dangerous hack.
        % Fundamentally, the problem is that this test is so sensitive.

        % 1.1102230246251565404236316680908203125e-16

        ((Vp_Dot_Vcp > Vs_Dot_Vcp) or (abs(Vp_Dot_Vcp - Vs_Dot_Vcp) < 2.0e-16))
    end.


% Returns true if the point PointGP is on the great circle plane.
% Forms the scalar triple product of A,B,p and if the volume of the
% resulting parallelepiped is near zero the point PointGP is on the
% great circle plane.
edge_point_on_plane({PointFromGP, PointToGP}, PointGP) ->

    % Normal to the plane defined by PointFromGP and PointToGP
    Normal = normalize(robust_cross_product(PointFromGP, PointToGP)),
    Point = geog2cart(PointGP),

    % We expect the dot product of with normal with any vector in the plane to be zero
    W = dot_product(Normal, Point),

    (W =< ?TOLERANCE).


% Tests if Inner is within Outer box
within({IW, IS, IE, IN}, {OW, OS, OE, ON}) ->
    (IW >= OW) and (IS >= OS) and (IE =< OE) and (IN =< ON);
% Tests if a point is within a box
within({X, Y}, {OW, OS, OE, ON}) ->
    (X >= OW) and (Y >= OS) and (X =< OE) and (Y =< ON).


% Modulo function for negative and double values according to the following specification:
%
% "These functions shall return the floating-point remainder r= x- ny when y is non-zero. 
% The value n is the integral value nearest the exact value x/ y. When |n-x/y|= 1/2, 
% the value n is chosen to be even."
% http://pubs.opengroup.org/onlinepubs/009604599/functions/remainder.html
remainder(A, B) ->
    A_div_B = A / B,
    N = round(A_div_B),

    case (abs(N - A_div_B) == 0.5) of
    true ->
        A_div_B_Trunc = trunc(A_div_B),

        New_N = case ((abs(A_div_B_Trunc) rem 2) == 0) of
        true -> A_div_B_Trunc;
        false ->
            case (A_div_B >= 0) of
            true -> A_div_B_Trunc + 1;
            false -> A_div_B_Trunc - 1
            end
        end,
        A - New_N * B;
    false ->
        A - N * B
    end.

