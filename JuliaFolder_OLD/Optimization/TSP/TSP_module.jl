module TSP_module

export greedy_shortest_path

#function greedy_shortest_path(D,start_node::Int64)
function greedy_shortest_path(D,start_node)
    num_of_cities = size(D)[1]
    current_node = start_node
    cities_visited = fill(0,num_of_cities)
    cities_visited[1] = start_node
    #i = 2
    for i=2:num_of_cities
        # Find the cities we can go to
        all_cities = collect(1:num_of_cities)
        #deleteat!(all_cities,find(l->l!=0,cities_visited))
        deleteat!(all_cities,sort!(cities_visited[find(l->l!=0,cities_visited)]))
        println(all_cities)
        # Find the city which is the smallest distance from current one
        #D[current_node,all_cities]
        println("does this show??")
        cities_visited[i] = find(l->l==minimum(D[current_node,all_cities]),D[current_node,:])[1]
        println(find(l->l==minimum(D[current_node,all_cities]),D[current_node,:])[1])
        println(cities_visited)
        current_node = cities_visited[i]
        println(current_node)
    end

    return cities_visited

end


end
