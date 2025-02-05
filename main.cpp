//Projekt 2 Projektowanie efektwnych algorytmów
//Autor Cyprian Kozubek nr.indeksu 272959
//Problem komiwojazera: Metoda podziału i ograniczeń (B&B – branch and bound).
// Należy zastosować trzy metody przeszukiwania przestrzeni rozwiązań. Wybrano
//• przeszukiwanie wszerz (breadth first search),
//• przeszukiwanie wgłab (depth first search),
//• przeszukiwanie pod kątem najniższego kosztu(lowest cost / best first search)

#include <iostream>
#include <windows.h>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <random>

using namespace std;

//Elementy zajmujące się mierzeniem czasu
double t1 = 0.0;
double t2 = 0.0;
double t3 = 0.0;
double PCFreq = 0.0;
__int64 CounterStart = 0;
void StartCounter(){
    LARGE_INTEGER li;
    if( !QueryPerformanceFrequency( & li ) )
        cout << "QueryPerformanceFrequency failed!\n";
    PCFreq = double( li.QuadPart ) / 1000.0;
    QueryPerformanceCounter( & li );
    CounterStart = li.QuadPart;
}
double GetCounter(){
    LARGE_INTEGER li;
    QueryPerformanceCounter( & li );
    return double( li.QuadPart - CounterStart ) / PCFreq;
}
int ** Matrix = nullptr; //Tablica dwuwymiarowa dla grafu
int matrix_size; //Rozmiar grafu
// Funkcja alokująca pamięć dla macierzy
void allocateMatrix(int rows, int cols) {
    Matrix = new int*[rows];
    for (int i = 0; i < rows; i++) {
        Matrix[i] = new int[cols];
    }
}
//Wczytywanie danych z pliku
void loadGraphFromFile(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Unable to open file: " << filename << endl;
        exit(1);
    }
    int n;
    file >> n;
    matrix_size = n;
    allocateMatrix(n,n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
               file >> Matrix[i][j];
        }
    }
    file.close();
}
// Generowanie losowego grafu (asymetrycznego)
void generate_graph(int graphSize) {
    matrix_size = graphSize;
    allocateMatrix(graphSize, graphSize);
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> distrib(1, 100);
    for (int i = 0; i < graphSize; i++) {
        for (int j = 0; j < graphSize; j++) {
            if (i == j) {
                Matrix[i][j] = -1;
            } else {
                int liczba = distrib(gen);
                Matrix[i][j] = liczba;
            }
        }
    }
}
// Generowanie losowego grafu (symetrycznego)
void generate_graph_sym(int graphSize) {
    matrix_size = graphSize;
    allocateMatrix(graphSize, graphSize);
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> distrib(1, 100);
    for (int i = 0; i < graphSize; i++) {
        for (int j = 0; j < graphSize; j++) {
            if (i == j) {
                Matrix[i][j] = -1;
            } else if (i < j) {
                int liczba = distrib(gen);
                Matrix[i][j] = liczba;
                Matrix[j][i] = liczba;
            }
        }
    }
}
//Wyswietlanie macierzy
void print_matrix(int size){
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            cout<<Matrix[i][j]<<" ";
        }
        cout<<endl;
    }
}
//Struktury dla algorytmow
//Struktura dla elementu trasy w BFS i DFS
struct Node {
    int vertex;        // Bieżący wierzchołek
    int* path;         // Aktualna trasa
    int cost;          // Koszt dotychczasowej trasy
    int pathLength;    // Liczba odwiedzonych miast
    bool* visited;     // Tablica miast odwiedzonych
    Node(int v, int n, int c, int pLength) : vertex(v), cost(c), pathLength(pLength) {
        path = new int[n];
        visited = new bool[n]();
    }
    ~Node() {
        delete[] path;
        delete[] visited;
    }
};
Node* allocate_node(int vertex, int cost, int path_length) {
    Node* newNode = new Node(vertex, matrix_size, cost, path_length);
    fill_n(newNode->path, matrix_size, -1);
    fill_n(newNode->visited, matrix_size, false);
    return newNode;
}
// Struktura stosu
struct Stack {
    int top;
    int capacity;
    Node** array;
    Stack() : top(-1), capacity(1) {
        array = new Node*[capacity];
    }
    bool isEmpty() {
        return top == -1;
    }
    void push(Node* item) {
        if (top + 1 == capacity) {
            capacity *= 2;
            Node** newArray = new Node*[capacity];
            copy(array, array + top + 1, newArray);
            delete[] array;
            array = newArray;
        }
        array[++top] = item;
    }
    Node* pop() {
        if (isEmpty()) return nullptr;
        return array[top--];
    }
    ~Stack() {
        for (int i = 0; i <= top; ++i) {
            delete array[i];
        }
        delete[] array;
    }
};
//Struktura kolejki priorytetowej
struct PriorityQueue {
    int size;
    int capacity;
    Node** array;
    PriorityQueue() : size(0), capacity(1) {
        array = new Node*[capacity];
    }
    bool isEmpty() {
        return size == 0;
    }
    void enqueue(Node* item) {
        if (size == capacity) {
            capacity *= 2;
            Node** newArray = new Node*[capacity];
            copy(array, array + size, newArray);
            delete[] array;
            array = newArray;
        }
        int i = size - 1;
        while (i >= 0 && array[i]->cost > item->cost) {
            array[i + 1] = array[i];
            i--;
        }
        array[i + 1] = item;
        size++;
    }
    Node* dequeue() {
        if (isEmpty()) return nullptr;
        Node* item = array[0];
        for (int i = 1; i < size; ++i) {
            array[i - 1] = array[i];
        }
        size--;
        return item;
    }
    // Przycinanie kolejki do określonego rozmiaru
    void trimToSize(int maxSize) {
        if (size > maxSize) {
            for (int i = maxSize; i < size; ++i) {
                delete array[i];
            }
            size = maxSize;
        }
    }
    ~PriorityQueue() {
        for (int i = 0; i < size; ++i) {
            delete array[i];
        }
        delete[] array;
    }
};
//Wyliczanie dolnego ograniczenia
int calculateLowerBound(int** distanceMatrix, int numVertices, Node* node) {
    int lowerBound = node->cost;
    for (int i = 0; i < numVertices; ++i) {
        if (!node->visited[i]) {
            int min1 = INT_MAX, min2 = INT_MAX;
            for (int j = 0; j < numVertices; ++j) {
                if (i != j && distanceMatrix[i][j] > 0) {
                    if (distanceMatrix[i][j] < min1) {
                        min2 = min1;
                        min1 = distanceMatrix[i][j];
                    } else if (distanceMatrix[i][j] < min2) {
                        min2 = distanceMatrix[i][j];
                    }
                }
            }
            if (min1 != INT_MAX) lowerBound += min1;
            if (min2 != INT_MAX) lowerBound += min2;
        }
    }
    return lowerBound / 2;
}
//Wyliczanie gornego ograniczenia
int calculateUpperBound(int** distanceMatrix, int numVertices, int currentBest) {
    int upperBound = currentBest;
    return upperBound;
}
//*************************Funkcje dla asymetrycznych grafow*****************************************************************
//BFS - Breadth First Search, przeszukiwanie wszerz
bool bfs(int** distanceMatrix, int numVertices, int source, int* best_route, int& best_cost) {
    PriorityQueue queue;
    Node* startNode = allocate_node(source, 0, 1);
    startNode->path[0] = source;
    startNode->visited[source] = true;
    queue.enqueue(startNode);
    int beamWidth = matrix_size*100;
    while (!queue.isEmpty()) {
        if (queue.size > beamWidth) {
            queue.trimToSize(beamWidth);
        }
        Node* currentNode = queue.dequeue();
        if (!currentNode) continue;
        if (currentNode->pathLength == numVertices) {
            int returnCost = distanceMatrix[currentNode->vertex][source];
            if (returnCost > 0) {
                int totalCost = currentNode->cost + returnCost;
                if (totalCost < best_cost) {
                    best_cost = totalCost;
                    copy(currentNode->path, currentNode->path + numVertices, best_route);
                    best_route[numVertices] = source;
                }
            }
            delete currentNode;
            continue;
        }
        for (int v = 0; v < numVertices; ++v) {
            if (!currentNode->visited[v] && distanceMatrix[currentNode->vertex][v] > 0) {
                int newCost = currentNode->cost + distanceMatrix[currentNode->vertex][v];
                Node* childNode = allocate_node(v, newCost, currentNode->pathLength + 1);
                copy(currentNode->path, currentNode->path + currentNode->pathLength, childNode->path);
                childNode->path[currentNode->pathLength] = v;
                copy(currentNode->visited, currentNode->visited + numVertices, childNode->visited);
                childNode->visited[v] = true;
                int lowerBound = calculateLowerBound(distanceMatrix, numVertices, childNode);
                if (newCost + lowerBound < best_cost) {
                    queue.enqueue(childNode);
                } else {
                    delete childNode;
                }
            }
        }
        delete currentNode;
    }
    return best_cost != INT_MAX;
}
//Algorytm Branch and Bound - BFS dla TSP
void branch_bound_BFS() {
    int best_cost = INT_MAX;
    int* best_route = new int[matrix_size + 1];
    int upperBound = calculateUpperBound(Matrix, matrix_size, best_cost);
    for (int start_vertex = 0; start_vertex < matrix_size; ++start_vertex) {
        int current_cost = upperBound;
        int* current_route = new int[matrix_size + 1];
        if (bfs(Matrix, matrix_size, start_vertex, current_route, current_cost)) {
            if (current_cost < best_cost) {
                best_cost = current_cost;
                copy(current_route, current_route + matrix_size + 1, best_route);
            }
        }
        delete[] current_route;
    }
    // Zapis wyników
    ofstream outfile("wyniki_branch_bound_BFS.txt");
    outfile << "Najmniejsza droga: " << best_cost << "\nTrasa: ";
    for (int i = 0; i < matrix_size; ++i) {
        outfile << best_route[i] << " -> ";
    }
    outfile << best_route[0] << "\n";
    outfile.close();
    // Wyświetlenie wyników
    cout << "Najmniejsza droga (B&B BFS): " << best_cost << "\nTrasa: ";
    for (int i = 0; i < matrix_size; ++i) {
        cout << best_route[i] << " -> ";
    }
    cout << best_route[0] << endl;
    delete[] best_route;
}
//DFS - Depth First Search, przeszukiwanie w glab
bool dfs(int** distanceMatrix, int numVertices, int source, int* best_route, int& best_cost) {
    Stack stack;
    Node* startNode = allocate_node(source, 0, 1);
    startNode->path[0] = source;
    startNode->visited[source] = true;
    stack.push(startNode);
    best_cost = INT_MAX;
    while (!stack.isEmpty()) {
        Node* currentNode = stack.pop();
        if (!currentNode) continue;
        if (currentNode->pathLength == numVertices) {
            int returnCost = distanceMatrix[currentNode->vertex][source];
            if (returnCost > 0) {
                int totalCost = currentNode->cost + returnCost;
                if (totalCost < best_cost) {
                    best_cost = totalCost;
                    copy(currentNode->path, currentNode->path + numVertices, best_route);
                    best_route[numVertices] = source;
                }
            }
            delete currentNode;
            continue;
        }
        for (int v = 0; v < numVertices; ++v) {
            if (!currentNode->visited[v] && distanceMatrix[currentNode->vertex][v] > 0) {
                int newCost = currentNode->cost + distanceMatrix[currentNode->vertex][v];
                Node* childNode = allocate_node(v, newCost, currentNode->pathLength + 1);
                copy(currentNode->path, currentNode->path + currentNode->pathLength, childNode->path);
                childNode->path[currentNode->pathLength] = v;
                copy(currentNode->visited, currentNode->visited + numVertices, childNode->visited);
                childNode->visited[v] = true;
                int lowerBound = calculateLowerBound(distanceMatrix, numVertices, childNode);
                if (newCost + lowerBound < best_cost) {
                    stack.push(childNode);
                } else {
                    delete childNode;
                }
            }
        }
        delete currentNode;
    }
    return best_cost != INT_MAX;
}
//Algorytm Branch and Bound - DFS dla TSP
void branch_bound_DFS() {
    int best_cost = INT_MAX;
    int* best_route = new int[matrix_size + 1];
    for (int start_vertex = 0; start_vertex < matrix_size; ++start_vertex) {
        int current_cost = INT_MAX;
        int* current_route = new int[matrix_size + 1];
        if (dfs(Matrix, matrix_size, start_vertex, current_route, current_cost)) {
            if (current_cost < best_cost) {
                best_cost = current_cost;
                copy(current_route, current_route + matrix_size + 1, best_route);
            }
        }
        delete[] current_route;
    }
    //Zapis do pliku
    ofstream outfile("wyniki_branch_bound_DFS.txt");
    outfile << "Najmniejsza droga: " << best_cost << "\nTrasa: ";
    for (int i = 0; i < matrix_size; ++i) {
        outfile << best_route[i] << " -> ";
    }
    outfile << best_route[0] << "\n";
    outfile.close();
    //Wyswietlanie wyniku
    cout << "Najmniejsza droga (B&B DFS): " << best_cost << "\nTrasa: ";
    for (int i = 0; i < matrix_size; ++i) {
        cout << best_route[i] << " -> ";
    }
    cout << best_route[0] << endl;
    delete[] best_route;
}
// Branch and Bound - Lowest Cost dla TSP
void branch_bound_lowest_cost() {
    int global_best_cost = INT_MAX;
    int* global_best_route = new int[matrix_size + 1];
    for (int start_vertex = 0; start_vertex < matrix_size; ++start_vertex) {
        PriorityQueue pq;
        int best_cost = global_best_cost;
        int* best_route = new int[matrix_size + 1];
        Node* startNode = allocate_node(start_vertex, 0, 1);
        startNode->path[0] = start_vertex;
        startNode->visited[start_vertex] = true;
        pq.enqueue(startNode);
        int beamWidth = matrix_size * 100;
        int iterationsWithoutImprovement = 0;
        while (!pq.isEmpty()) {
            pq.trimToSize(beamWidth);
            Node* currentNode = pq.dequeue();
            if (!currentNode) continue;
            if (currentNode->pathLength == matrix_size) {
                int returnCost = Matrix[currentNode->vertex][start_vertex];
                if (returnCost > 0) {
                    int totalCost = currentNode->cost + returnCost;
                    if (totalCost < best_cost) {
                        best_cost = totalCost;
                        copy(currentNode->path, currentNode->path + matrix_size, best_route);
                        best_route[matrix_size] = start_vertex;
                        iterationsWithoutImprovement = 0;
                    }
                }
                delete currentNode;
                continue;
            }
            for (int v = 0; v < matrix_size; ++v) {
                if (!currentNode->visited[v] && Matrix[currentNode->vertex][v] > 0) {
                    int newCost = currentNode->cost + Matrix[currentNode->vertex][v];
                    int lowerBound = calculateLowerBound(Matrix, matrix_size, currentNode);
                    if (newCost + lowerBound >= best_cost) continue;
                    Node* childNode = allocate_node(v, newCost, currentNode->pathLength + 1);
                    copy(currentNode->path, currentNode->path + currentNode->pathLength, childNode->path);
                    childNode->path[currentNode->pathLength] = v;
                    copy(currentNode->visited, currentNode->visited + matrix_size, childNode->visited);
                    childNode->visited[v] = true;
                    pq.enqueue(childNode);
                }
            }
            delete currentNode;
            // Dynamiczne dostosowanie beamWidth
            if (pq.size < beamWidth / 2) {
                beamWidth = min(matrix_size * 1000, beamWidth * 2);
            } else if (iterationsWithoutImprovement > 50) {
                beamWidth = max(matrix_size, beamWidth / 2);
                iterationsWithoutImprovement = 0;
            } else {
                iterationsWithoutImprovement++;
            }
        }
        if (best_cost < global_best_cost) {
            global_best_cost = best_cost;
            copy(best_route, best_route + matrix_size + 1, global_best_route);
        }
        delete[] best_route;
    }
    // Zapis wyników do pliku
    ofstream outfile("wyniki_branch_bound_lowest_cost.txt");
    outfile << "Najmniejsza droga: " << global_best_cost << "\nTrasa: ";
    for (int i = 0; i < matrix_size; ++i) {
        outfile << global_best_route[i] << " -> ";
    }
    outfile << global_best_route[0] << "\n";
    outfile.close();
    // Wyświetlenie wyników
    cout << "Najmniejsza droga (B&B Lowest Cost): " << global_best_cost << "\nTrasa: ";
    for (int i = 0; i < matrix_size; ++i) {
        cout << global_best_route[i] << " -> ";
    }
    cout << global_best_route[0] << endl;
    delete[] global_best_route;
}
//**************************Funkcje dla symetrycznych grafow*****************************************************************
bool bfs_sym(int** distanceMatrix, int numVertices, int source, int* best_route, int& best_cost) {
    PriorityQueue queue;
    Node* startNode = allocate_node(source, 0, 1);
    startNode->path[0] = source;
    startNode->visited[source] = true;
    queue.enqueue(startNode);
    int beamWidth = numVertices * 50;
    while (!queue.isEmpty()) {
        if (queue.size > beamWidth) {
            queue.trimToSize(beamWidth);
        }
        Node* currentNode = queue.dequeue();
        if (!currentNode) continue;
        if (currentNode->pathLength == numVertices) {
            int returnCost = distanceMatrix[currentNode->vertex][source];
            if (returnCost > 0) {
                int totalCost = currentNode->cost + returnCost;
                if (totalCost < best_cost) {
                    best_cost = totalCost;
                    copy(currentNode->path, currentNode->path + numVertices, best_route);
                    best_route[numVertices] = source;
                }
            }
            delete currentNode;
            continue;
        }
        for (int v = 0; v < numVertices; ++v) {
            if (!currentNode->visited[v] && v != currentNode->vertex && distanceMatrix[currentNode->vertex][v] != -1) {
                int newCost = currentNode->cost + distanceMatrix[currentNode->vertex][v];
                Node* childNode = allocate_node(v, newCost, currentNode->pathLength + 1);
                copy(currentNode->path, currentNode->path + currentNode->pathLength, childNode->path);
                childNode->path[currentNode->pathLength] = v;
                copy(currentNode->visited, currentNode->visited + numVertices, childNode->visited);
                childNode->visited[v] = true;
                int lowerBound = calculateLowerBound(distanceMatrix, numVertices, childNode);

                if (newCost + lowerBound < best_cost) {
                    queue.enqueue(childNode);
                } else {
                    delete childNode;
                }
            }
        }
        delete currentNode;
    }
    return best_cost != INT_MAX;
}
bool dfs_sym(int** distanceMatrix, int numVertices, int source, int* best_route, int& best_cost) {
    Stack stack;
    Node* startNode = allocate_node(source, 0, 1);
    startNode->path[0] = source;
    startNode->visited[source] = true;
    stack.push(startNode);
    best_cost = INT_MAX;
    while (!stack.isEmpty()) {
        Node* currentNode = stack.pop();
        if (!currentNode) continue;
        if (currentNode->pathLength == numVertices) {
            int returnCost = distanceMatrix[currentNode->vertex][source];
            if (returnCost > 0) {
                int totalCost = currentNode->cost + returnCost;
                if (totalCost < best_cost) {
                    best_cost = totalCost;
                    copy(currentNode->path, currentNode->path + numVertices, best_route);
                    best_route[numVertices] = source;
                }
            }
            delete currentNode;
            continue;
        }
        for (int v = 0; v < numVertices; ++v) {
            if (!currentNode->visited[v] && v != currentNode->vertex && distanceMatrix[currentNode->vertex][v] != -1) {
                int newCost = currentNode->cost + distanceMatrix[currentNode->vertex][v];
                int lowerBound = calculateLowerBound(distanceMatrix, numVertices, currentNode);
                if (newCost + lowerBound >= best_cost) {
                    continue;
                }
                Node* childNode = allocate_node(v, newCost, currentNode->pathLength + 1);
                copy(currentNode->path, currentNode->path + currentNode->pathLength, childNode->path);
                childNode->path[currentNode->pathLength] = v;
                copy(currentNode->visited, currentNode->visited + numVertices, childNode->visited);
                childNode->visited[v] = true;
                stack.push(childNode);
            }
        }
        delete currentNode;
    }
    return best_cost != INT_MAX;
}
void branch_bound_BFS_sym() {
    int best_cost = INT_MAX;
    int* best_route = new int[matrix_size + 1];
    for (int start_vertex = 0; start_vertex < matrix_size; ++start_vertex) {
        int current_cost = INT_MAX;
        int* current_route = new int[matrix_size + 1];
        if (bfs_sym(Matrix, matrix_size, start_vertex, current_route, current_cost)) {
            if (current_cost < best_cost) {
                best_cost = current_cost;
                copy(current_route, current_route + matrix_size + 1, best_route);
            }
        }
        delete[] current_route;
    }
    // Zapis wyników
    ofstream outfile("wyniki_branch_bound_BFS.txt");
    outfile << "Najmniejsza droga: " << best_cost << "\nTrasa: ";
    for (int i = 0; i < matrix_size; ++i) {
        outfile << best_route[i] << " -> ";
    }
    outfile << best_route[0] << "\n";
    outfile.close();
    cout << "Najmniejsza droga (B&B BFS Sym): " << best_cost << "\nTrasa: ";
    for (int i = 0; i < matrix_size; ++i) {
        cout << best_route[i] << " -> ";
    }
    cout << best_route[0] << endl;
    delete[] best_route;
}
void branch_bound_DFS_sym() {
    int best_cost = INT_MAX;
    int* best_route = new int[matrix_size + 1];
    for (int start_vertex = 0; start_vertex < matrix_size; ++start_vertex) {
        int current_cost = INT_MAX;
        int* current_route = new int[matrix_size + 1];
        if (dfs_sym(Matrix, matrix_size, start_vertex, current_route, current_cost)) {
            if (current_cost < best_cost) {
                best_cost = current_cost;
                copy(current_route, current_route + matrix_size + 1, best_route);
            }
        }
        delete[] current_route;
    }
    // Zapis wyników
    ofstream outfile("wyniki_branch_bound_DFS.txt");
    outfile << "Najmniejsza droga: " << best_cost << "\nTrasa: ";
    for (int i = 0; i < matrix_size; ++i) {
        outfile << best_route[i] << " -> ";
    }
    outfile << best_route[0] << "\n";
    outfile.close();
    cout << "Najmniejsza droga (B&B DFS Sym): " << best_cost << "\nTrasa: ";
    for (int i = 0; i < matrix_size; ++i) {
        cout << best_route[i] << " -> ";
    }
    cout << best_route[0] << endl;
    delete[] best_route;
}
void branch_bound_lowest_cost_sym() {
    int global_best_cost = INT_MAX;
    int* global_best_route = new int[matrix_size + 1];
    for (int start_vertex = 0; start_vertex < matrix_size; ++start_vertex) {
        PriorityQueue pq;
        int best_cost = global_best_cost;
        int* best_route = new int[matrix_size + 1];
        Node* startNode = allocate_node(start_vertex, 0, 1);
        startNode->path[0] = start_vertex;
        startNode->visited[start_vertex] = true;
        pq.enqueue(startNode);
        int beamWidth = matrix_size * 50;
        int iterationsWithoutImprovement = 0;
        while (!pq.isEmpty()) {
            pq.trimToSize(beamWidth);
            Node* currentNode = pq.dequeue();
            if (!currentNode) continue;
            if (currentNode->pathLength == matrix_size) {
                int returnCost = Matrix[start_vertex][currentNode->vertex];
                if (returnCost > 0) {
                    int totalCost = currentNode->cost + returnCost;
                    if (totalCost < best_cost) {
                        best_cost = totalCost;
                        copy(currentNode->path, currentNode->path + matrix_size, best_route);
                        best_route[matrix_size] = start_vertex;
                    }
                }
                delete currentNode;
                continue;
            }
            for (int v = 0; v < matrix_size; ++v) {
                if (v == currentNode->vertex || currentNode->visited[v] || Matrix[currentNode->vertex][v] == -1)
                    continue;
                int newCost = currentNode->cost + Matrix[currentNode->vertex][v];
                int lowerBound = calculateLowerBound(Matrix, matrix_size, currentNode);
                if (newCost + lowerBound >= best_cost) continue;
                Node* childNode = allocate_node(v, newCost, currentNode->pathLength + 1);
                copy(currentNode->path, currentNode->path + currentNode->pathLength, childNode->path);
                childNode->path[currentNode->pathLength] = v;
                copy(currentNode->visited, currentNode->visited + matrix_size, childNode->visited);
                childNode->visited[v] = true;
                pq.enqueue(childNode);
            }
            delete currentNode;
            //Dynamiczne ustawianie beamWidth
            if (pq.size < beamWidth / 2) {
                beamWidth = min(matrix_size * 500, beamWidth * 2);
            } else if (iterationsWithoutImprovement > 50) {
                beamWidth = max(matrix_size, beamWidth / 2);
                iterationsWithoutImprovement = 0;
            } else {
                iterationsWithoutImprovement++;
            }
        }
        if (best_cost < global_best_cost) {
            global_best_cost = best_cost;
            copy(best_route, best_route + matrix_size + 1, global_best_route);
        }
        delete[] best_route;
    }

    ofstream outfile("wyniki_branch_bound_lowest_cost_sym.txt");
    outfile << "Najmniejsza droga: " << global_best_cost << "\nTrasa: ";
    for (int i = 0; i < matrix_size; ++i) {
        outfile << global_best_route[i] << " -> ";
    }
    outfile << global_best_route[0] << "\n";
    outfile.close();
    cout << "Najmniejsza droga (B&B Lowest Cost Sym): " << global_best_cost << "\nTrasa: ";
    for (int i = 0; i < matrix_size; ++i) {
        cout << global_best_route[i] << " -> ";
    }
    cout << global_best_route[0] << endl;
    delete[] global_best_route;
}

// Funkcja do wczytania danych z pliku konfiguracyjnego
bool loadConfig(const string& configFile, int& graphType, int& algorithmType, int& load_graph, string& graphFile, int& repeatCount, int& graphSize, int& displayGraph) {
    ifstream config(configFile);
    if (!config.is_open()) {
        cout << "Blad otwarcia pliku konfiguracyjnego!" << endl;
        return false;
    }
    string line;
    while (getline(config, line)) {
        stringstream ss(line);
        string key, value;
        if (line.empty() || line[0] == '#') continue;
        getline(ss, key, '=');
        getline(ss, value);
        key.erase(0, key.find_first_not_of(" \t"));
        key.erase(key.find_last_not_of(" \t") + 1);
        value.erase(0, value.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t") + 1);
        if (key == "graph_type") {
            graphType = stoi(value);
        } else if (key == "algorithm") {
            algorithmType = stoi(value);
        } else if (key == "load_graph") {
            load_graph = stoi(value);
        } else if (key == "graph_file") {
            graphFile = value;
        } else if (key == "repeat_count") {
            repeatCount = stoi(value);
        } else if (key == "graph_size") {
            graphSize = stoi(value);
        } else if (key == "display_graph") {
            displayGraph = stoi(value);
        }
    }
    config.close();
    return true;
}
// Menu wyboru typu generowanego grafu
void Menu_Type(int graphType, int graphSize) {
    if (graphType == 1) {
        generate_graph(graphSize);
    } else if (graphType == 2) {
        generate_graph_sym(graphSize);
    }
}
// Menu wyboru algorytmów
void Menu_Alg(int algorithmType, int repeatCount, int graphType) {
        if (algorithmType == 0 || algorithmType == 1) {
            StartCounter();
            if (graphType == 1) {
                for (int i = 0; i < repeatCount; ++i) {
                    branch_bound_BFS();
                }
            } else if (graphType == 2) {
                for (int i = 0; i < repeatCount; ++i) {
                    branch_bound_BFS_sym();
                }
            }
            t1 = GetCounter();
            cout << "Czas dzialania BFS w ms: " << t1/repeatCount << endl;
        }
        if (algorithmType == 0 || algorithmType == 2) {
            StartCounter();
            if (graphType == 1) {
                for (int i = 0; i < repeatCount; ++i) {
                    branch_bound_DFS();
                }
            } else if (graphType == 2) {
                for (int i = 0; i < repeatCount; ++i) {
                    branch_bound_DFS_sym();
                }
            }
            t2 = GetCounter();
            cout << "Czas dzialania DFS w ms: " << t2/repeatCount << endl;
        }
        if (algorithmType == 0 || algorithmType == 3) {
            StartCounter();
            if (graphType == 1) {
                for (int i = 0; i < repeatCount; ++i) {
                    branch_bound_lowest_cost();
                }
            } else if (graphType == 2) {
                for (int i = 0; i < repeatCount; ++i) {
                    branch_bound_lowest_cost_sym();
                }
            }
            t3 = GetCounter();
            cout << "Czas dzialania Lowest Cost w ms: " << t3/repeatCount << endl;
        }
    }
// Menu główne programu
void Menu() {
    int graphType = 0, algorithmType = 0, load_graph = 0, repeatCount = 1, graphSize = 0, displayGraph = 0;
    string graphFile;
    //Należy podać ścieżkę bezwzględną pliku
    if (!loadConfig(R"()",
                    graphType, algorithmType, load_graph, graphFile, repeatCount, graphSize, displayGraph)) {
        return;
    }
    if (load_graph == 1 && !graphFile.empty()) {
        loadGraphFromFile(graphFile);
    } else if (load_graph == 2) {
        Menu_Type(graphType, graphSize);
    } else {
        cout << "Nieprawidlowy wybor opcji wczytywania grafu!" << endl;
        return;
    }
    if (displayGraph == 1) {
        cout << "Aktualny graf:" << endl;
        print_matrix(matrix_size);
    }
    Menu_Alg(algorithmType, repeatCount, graphType);
}
int main() {
    Menu();
    return 0;
}
