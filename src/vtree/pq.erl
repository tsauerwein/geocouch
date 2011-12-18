-module(pq).

% non-optimal priority queue which uses the module gb_trees

-export([empty/0, add/3, isEmpty/1, takeMin/1]).

empty() ->
    gb_trees:empty().

add(Key, Val, Tree) ->
    Result = gb_trees:lookup(Key, Tree),
    
    case Result of
    none -> 
        gb_trees:insert(Key, [Val], Tree);
    {value, OldVals} ->
        gb_trees:enter(Key, [Val | OldVals], Tree)
    end.

isEmpty(Tree) ->
    gb_trees:is_empty(Tree).
    
takeMin(Tree) ->
    {Key, Vals, NewTree} = gb_trees:take_smallest(Tree),
    
    case Vals of
    [Val | []] -> 
        {Val, NewTree};
    [Val | RemainingVals] ->
        {Val, gb_trees:enter(Key, RemainingVals, Tree)}
    end.
