#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <sstream>
#include <utility>
#include <cmath>
#include <algorithm>

class k_heap {
private:
  int* heap_arr;
  int* positions; 
  std::vector<long long> dist;
  int capacity;
  int current_size;
  int k;

  void swap(int i, int j) {
    int temp = heap_arr[i];
    heap_arr[i] = heap_arr[j];
    heap_arr[j] = temp;

    positions[heap_arr[i]] = i;
    positions[heap_arr[j]] = j;
  }

  int parent(int i) const {
    return (i - 1) / k;
  }
  
  int child(int i, int j) const {
    return k * i + j + 1;
  }

  bool is_leaf(int i) {
    return child(i, 0) >= current_size;
  }

  void heapify_up(int i) {
    if (root(i)) return;

    if (key(parent(i)) > key(i)) {
      swap(i, parent(i));
      heapify_up(parent(i));
    }
  }

  void heapify_down(int i) {
    if (is_leaf(i)) return;

    int smallest_child = child(i, 0);
    for (int j = 1; j < k; j++) {
      int child_idx = child(i, j);
      if (child_idx < current_size && key(child_idx) < key(smallest_child)) {
        smallest_child = child_idx;
      }
    }

    if (key(smallest_child) < key(i)) {
      swap(i, smallest_child);
      heapify_down(smallest_child);
    }
  }

public:
  k_heap(int cap, int k_value) {
    if (k_value < 2) {
      throw std::invalid_argument("k must be at least 2");
    }
    capacity = cap;
    heap_arr = new int[capacity];
    positions = new int[capacity];
    dist.resize(capacity, std::numeric_limits<long long>::max()); 
    current_size = 0;
    k = k_value; 

    for (int i = 0; i < capacity; i++) {
      positions[i] = -1;
    }
  }

  ~k_heap() {
    delete[] heap_arr;
    delete[] positions;
  }

  long long key(int i) {  
    return dist[heap_arr[i]];
  }

  bool root(int i) {
    return i == 0;
  }

  void insert(int vertex) {
    heap_arr[current_size] = vertex;
    positions[vertex] = current_size;
    heapify_up(current_size);
    current_size++;
  }

  void deleteKey(int i) {
    if (i >= current_size) {
      throw std::out_of_range("Index out of range");
    }

    int vertex_to_remove = heap_arr[i];
    heap_arr[i] = heap_arr[current_size - 1];
    positions[heap_arr[i]] = i;
    positions[vertex_to_remove] = -1;
    current_size--;

    if (i > 0 && key(i) < key(parent(i))) {
      heapify_up(i);
    } else {
      heapify_down(i);
    }
  }

  void update(int vertex, long long new_value) {  
    if (positions[vertex] == -1) {
      dist[vertex] = new_value;
      insert(vertex);
    } else {
      long long old_value = dist[vertex];
      dist[vertex] = new_value;
      
      int i = positions[vertex];
      if (new_value < old_value) {
        heapify_up(i);
      } else if (new_value > old_value) {
        heapify_down(i);
      }
    }
  }

  void deleteMin() {
    if (current_size <= 0) {
      throw std::runtime_error("Heap is empty");
    }
    deleteKey(0);
  }

  int getMin() {
    if (current_size <= 0) {
      throw std::runtime_error("Heap is empty");
    }
    return heap_arr[0];
  }

  bool isEmpty() {
    return current_size == 0;
  }

  int size() {
    return current_size;
  }
  
  long long getDistance(int vertex) {  
    return dist[vertex];
  }
  
  void setDistance(int vertex, long long distance) {  
    dist[vertex] = distance;
  }
  
  bool isInHeap(int vertex) {
    return positions[vertex] != -1;
  }
  
  void clear() {
    current_size = 0;
    std::fill(dist.begin(), dist.end(), std::numeric_limits<long long>::max());
    std::fill(positions, positions + capacity, -1);
  }
};

typedef std::vector<std::vector<std::pair<int, long long>>> graph;  

void dijkstra(const graph& g, const std::vector<int>& sources,  
                       std::vector<long long>& distances,                
                       std::vector<int>& predecessors,                    
                       int k_value) {
  int n = g.size() - 1;
  
  k_heap Q(n + 1, k_value);
  
  for (int v = 1; v <= n; v++) {
    Q.setDistance(v, std::numeric_limits<long long>::max());
  }
  
  for (int src : sources) {
    Q.setDistance(src, 0);
    Q.insert(src);
  }
  
  std::vector<bool> visited(n + 1, false);
  predecessors.assign(n + 1, -1);  
  
  while (!Q.isEmpty()) {
    int v = Q.getMin();
    Q.deleteMin();
    
    visited[v] = true;
        
    for (const auto& edge : g[v]) {
      int u = edge.first;
      long long dvu = edge.second; 
      
      if (!visited[u]) {
        long long dv = Q.getDistance(v);
        long long du = Q.getDistance(u);
        
        if (du == std::numeric_limits<long long>::max()) {
          Q.setDistance(u, dv + dvu);
          Q.insert(u);
          predecessors[u] = v;  
        }
        else if (dv + dvu < du) {
          Q.update(u, dv + dvu);
          predecessors[u] = v; 
        }
      }
    }
  }
  
  for (int v = 1; v <= n; v++) {
    distances[v] = Q.getDistance(v);
  }
}

class Hungarian {
private:
  int n;
  std::vector<std::vector<long long>> weights;
  std::vector<long long> potentials;   
  std::vector<int> matching;           
  
  graph residual_graph;
  std::vector<long long> distances;
  std::vector<int> predecessors;
  std::vector<int> free_S, free_T;

  void init_potentials() {
    long long W = 0;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        W = std::max(W, weights[i][j]);
      }
    }
    
    potentials.resize(2 * n + 1);
    for (int v = 1; v <= n; v++) {
      potentials[v] = 0;      
    }
    for (int v = n + 1; v <= 2 * n; v++) {
      potentials[v] = W;        
    }
  }

  void find_free_vertices() {
    free_S.clear();
    free_T.clear();
    
    for (int v = 1; v <= n; v++) {
      if (matching[v] == -1) {
        free_S.push_back(v);
      }
    }
    
    for (int v = n + 1; v <= 2 * n; v++) {
      if (matching[v] == -1) {
        free_T.push_back(v);
      }
    }
  }

  void build_residual_graph() {
    residual_graph.assign(2 * n + 1, std::vector<std::pair<int, long long>>());
    
    for (int i = 1; i <= n; i++) {
      for (int j = 1; j <= n; j++) {
        int v = j + n;
        long long transformed_weight = weights[i-1][j-1] - potentials[v] + potentials[i];
        residual_graph[i].push_back(std::make_pair(v, transformed_weight));
      }
    }
    
    for (int u = n + 1; u <= 2 * n; u++) {
      if (matching[u] != -1) {
        int v = matching[u];
        long long reverse_weight = -weights[v-1][u-n-1] - potentials[v] + potentials[u];
        residual_graph[u].push_back(std::make_pair(v, reverse_weight));
      }
    }
  }

  void get_path(int target, std::vector<int>& path) {
    path.clear();
    int current = target;
    while (current != -1) {
      path.push_back(current);
      current = predecessors[current];
    }
    std::reverse(path.begin(), path.end());
  }

  void apply_path(const std::vector<int>& path) {
    for (int i = 0; i < path.size() - 1; i += 2) {
      int u = path[i];
      int v = path[i + 1];
      matching[u] = v;
      matching[v] = u;
    }
  }

  void update_potentials() {
    for (int v = 1; v <= 2 * n; v++) {
      if (distances[v] < std::numeric_limits<long long>::max()) {
        potentials[v] += distances[v];
      }
    }
  }

public:
  Hungarian(int size, const std::vector<std::vector<long long>>& w) 
    : n(size), weights(w), matching(2 * size + 1, -1) {
    init_potentials();
    distances.resize(2 * n + 1);
    predecessors.resize(2 * n + 1);
  }

  long long match() {
    int matched = 0;
    
    while (matched < n) {
      find_free_vertices();
      
      if (free_S.empty()) break;
      
      build_residual_graph();
      
      dijkstra(residual_graph, free_S, distances, predecessors, 2);
      
      long long min_dist = std::numeric_limits<long long>::max();
      int target = -1;
      
      for (int t : free_T) {
        if (distances[t] < min_dist) {
          min_dist = distances[t];
          target = t;
        }
      }
      
      if (target == -1 || distances[target] == std::numeric_limits<long long>::max()) {
        break;
      }
      
      std::vector<int> path;
      get_path(target, path);
      apply_path(path);
      
      update_potentials();
      
      matched++;
    }
    
    long long total_cost = 0;
    for (int i = 1; i <= n; i++) {
      if (matching[i] != -1) {
        int j = matching[i] - n;
        total_cost += weights[i-1][j-1];
      }
    }
    
    return total_cost;
  }
};

int main() {
  int n;
  std::cin >> n;
  
  std::vector<std::vector<long long>> weights(n, std::vector<long long>(n));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      std::cin >> weights[i][j];
    }
  }
  
  Hungarian hungarian(n, weights);
  long long result = hungarian.match();
  
  std::cout << result << std::endl;
  
  return 0;
}
