# Dijkstra_algorithm
# Abstract
Dijkstra's algorithm is an algorithm for finding the shortest paths between nodes in a graph, which may represent, for example, road networks. For a given source node in the graph, the algorithm finds the shortest path between that node and every other. It can also be used for finding the shortest paths from a single node to a single destination node by stopping the algorithm once the shortest path to the destination node has been determined. 


In Dijkstra’s algorithm, two sets are maintained, one set contains list of vertices already included in SPT (Shortest Path Tree), other set contains vertices not yet included. The idea is to traverse all vertices of graph using BFS and use a Min Heap to store the vertices not yet included in SPT (or the vertices for which shortest distance is not finalized yet).  Min Heap is used as a priority queue to get the minimum distance vertex from set of not yet included vertices.
# Problem Statement
Dijkstra’s algorithm: Given a grid of dimensions n × n. Every corner of a square is a junction (gate) which can be either closed or open. The goal is to reach from a specified starting position to a specified ending position by taking the shortest path (if exists). At every junction, a decision has to be made to take a turn (left/right) or continue moving straight.

Expected Output: Sequence required to traverse the grid consisting of right (R), left(L), Straight (S)
