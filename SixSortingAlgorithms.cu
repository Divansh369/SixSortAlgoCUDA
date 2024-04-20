#include <iostream>
#include <chrono>
#include <vector>
#include <random>
#include <cuda_runtime.h>
#include <iomanip>

using namespace std;

// Function to generate a random array
vector<int> generate_random_array(int size) {
    vector<int> arr(size);
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dis(-10000 , 10000);
    for (int i = 0; i < size; ++i) {
        arr[i] = dis(gen);
    }
    return arr;
}

// Quick Sort Serial
void quick_sort_serial(vector<int>& arr, int low, int high) {
    if (low < high) {
        int pivot = arr[high];
        int i = low - 1;
        for (int j = low; j <= high - 1; j++) {
            if (arr[j] < pivot) {
                i++;
                swap(arr[i], arr[j]);
            }
        }
        swap(arr[i + 1], arr[high]);
        int pi = i + 1;

        quick_sort_serial(arr, low, pi - 1);
        quick_sort_serial(arr, pi + 1, high);
    }
}

// Merge Sort Serial
void merge(vector<int>& arr, int l, int m, int r) {
    int n1 = m - l + 1;
    int n2 = r - m;

    vector<int> L(n1), R(n2);
    for (int i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (int j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];

    int i = 0;
    int j = 0;
    int k = l;
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            i++;
        }
        else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }

    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
}

void merge_sort_serial(vector<int>& arr, int l, int r) {
    if (l >= r) {
        return;
    }
    int m = l + (r - l) / 2;
    merge_sort_serial(arr, l, m);
    merge_sort_serial(arr, m + 1, r);
    merge(arr, l, m, r);
}

// Heap Sort Serial
void heapify(vector<int>& arr, int n, int i) {
    int largest = i;
    int l = 2 * i + 1;
    int r = 2 * i + 2;

    if (l < n && arr[l] > arr[largest])
        largest = l;

    if (r < n && arr[r] > arr[largest])
        largest = r;

    if (largest != i) {
        swap(arr[i], arr[largest]);
        heapify(arr, n, largest);
    }
}

void heap_sort_serial(vector<int>& arr) {
    int n = arr.size();

    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);

    for (int i = n - 1; i > 0; i--) {
        swap(arr[0], arr[i]);
        heapify(arr, i, 0);
    }
}

// Bubble Sort Serial
void bubble_sort_serial(vector<int>& arr) {
    int n = arr.size();
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < n - i - 1; j++) {
            if (arr[j] > arr[j + 1]) {
                swap(arr[j], arr[j + 1]);
            }
        }
    }
}

// Selection Sort Serial
void selection_sort_serial(vector<int>& arr) {
    int n = arr.size();
    for (int i = 0; i < n - 1; i++) {
        int min_idx = i;
        for (int j = i + 1; j < n; j++) {
            if (arr[j] < arr[min_idx]) {
                min_idx = j;
            }
        }
        swap(arr[i], arr[min_idx]);
    }
}

// Insertion Sort Serial
void insertion_sort_serial(vector<int>& arr) {
    int n = arr.size();
    for (int i = 1; i < n; i++) {
        int key = arr[i];
        int j = i - 1;
        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}

// Bubble Sort Parallel (CUDA)
__global__ void bubble_sort_kernel(int* arr, int n) {
    int idx = threadIdx.x;
    if (idx < n) {
        for (int i = 0; i < n - 1; i++) {
            if (idx < n - i - 1 && arr[idx] > arr[idx + 1]) {
                int temp = arr[idx];
                arr[idx] = arr[idx + 1];
                arr[idx + 1] = temp;
            }
        }
    }
}

// Selection Sort Parallel (CUDA)
__global__ void selection_sort_kernel(int* arr, int n) {
    int idx = threadIdx.x;
    if (idx < n - 1) {
        int min_idx = idx;
        for (int j = idx + 1; j < n; j++) {
            if (arr[j] < arr[min_idx]) {
                min_idx = j;
            }
        }
        int temp = arr[idx];
        arr[idx] = arr[min_idx];
        arr[min_idx] = temp;
    }
}

// Insertion Sort Parallel (CUDA)
__global__ void insertion_sort_kernel(int* arr, int n) {
    int idx = threadIdx.x;
    if (idx > 0 && idx < n) {
        int key = arr[idx];
        int j = idx - 1;
        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}

// CUDA kernel for parallel Quick Sort
__device__ void quick_sort_kernel(int* arr, int low, int high) {
    if (low < high) {
        int pivot = arr[high];
        int i = low - 1;
        for (int j = low; j <= high - 1; j++) {
            if (arr[j] < pivot) {
                i++;
                int temp = arr[i];
                arr[i] = arr[j];
                arr[j] = temp;
            }
        }
        int temp = arr[i + 1];
        arr[i + 1] = arr[high];
        arr[high] = temp;

        int pi = i + 1;

        quick_sort_kernel(arr, low, pi - 1);
        quick_sort_kernel(arr, pi + 1, high);
    }
}

__global__ void parallel_quick_sort_kernel(int* arr, int size) {
    quick_sort_kernel(arr, 0, size - 1);
}

// CUDA kernel for parallel Merge Sort
__device__ void merge_kernel(int* arr, int l, int m, int r) {
    int n1 = m - l + 1;
    int n2 = r - m;

    int* L = new int[n1];
    int* R = new int[n2];

    for (int i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (int j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];

    int i = 0;
    int j = 0;
    int k = l;
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            i++;
        } else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }

    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }

    delete[] L;
    delete[] R;
}

__device__ void merge_sort_kernel(int* arr, int l, int r) {
    if (l < r) {
        int m = l + (r - l) / 2;
        merge_sort_kernel(arr, l, m);
        merge_sort_kernel(arr, m + 1, r);
        merge_kernel(arr, l, m, r);
    }
}

__global__ void parallel_merge_sort_kernel(int* arr, int size) {
    merge_sort_kernel(arr, 0, size - 1);
}

// CUDA kernel for parallel Heap Sort
__device__ void max_heapify_kernel(int* arr, int n, int i) {
    int largest = i;
    int l = 2 * i + 1;
    int r = 2 * i + 2;

    if (l < n && arr[l] > arr[largest])
        largest = l;

    if (r < n && arr[r] > arr[largest])
        largest = r;

    if (largest != i) {
        int temp = arr[i];
        arr[i] = arr[largest];
        arr[largest] = temp;

        max_heapify_kernel(arr, n, largest);
    }
}

__global__ void parallel_heap_sort_kernel(int* arr, int size) {
    for (int i = size / 2 - 1; i >= 0; i--)
        max_heapify_kernel(arr, size, i);

    for (int i = size - 1; i >= 0; i--) {
        int temp = arr[0];
        arr[0] = arr[i];
        arr[i] = temp;

        max_heapify_kernel(arr, i, 0);
    }
}

// Serial Bubble Sort Execution
void serial_bubble_sort_execution(const vector<int>& arr) {
    vector<int> serial_arr(arr);
    auto start_serial = chrono::steady_clock::now();
    bubble_sort_serial(serial_arr);
    auto end_serial = chrono::steady_clock::now();
    chrono::duration<double> elapsed_serial = end_serial - start_serial;
    // cout << "Serial Bubble Sort Execution Time: " << elapsed_serial.count() * 1000 << " milliseconds" << endl;
}

// Serial Selection Sort Execution
void serial_selection_sort_execution(const vector<int>& arr) {
    vector<int> serial_arr(arr);
    auto start_serial = chrono::steady_clock::now();
    selection_sort_serial(serial_arr);
    auto end_serial = chrono::steady_clock::now();
    chrono::duration<double> elapsed_serial = end_serial - start_serial;
    // cout << "Serial Selection Sort Execution Time: " << elapsed_serial.count() * 1000 << " milliseconds" << endl;
}

// Serial Insertion Sort Execution
void serial_insertion_sort_execution(const vector<int>& arr) {
    vector<int> serial_arr(arr);
    auto start_serial = chrono::steady_clock::now();
    insertion_sort_serial(serial_arr);
    auto end_serial = chrono::steady_clock::now();
    chrono::duration<double> elapsed_serial = end_serial - start_serial;
    // cout << "Serial Insertion Sort Execution Time: " << elapsed_serial.count() * 1000 << " milliseconds" << endl;
}

// Serial Quick Sort Execution
void serial_quick_sort_execution(const vector<int>& arr) {
    vector<int> serial_arr(arr);
    auto start_serial = chrono::steady_clock::now();
    quick_sort_serial(serial_arr, 0, serial_arr.size() - 1);
    auto end_serial = chrono::steady_clock::now();
    chrono::duration<double> elapsed_serial = end_serial - start_serial;
    // cout << "Serial Quick Sort Execution Time: " << elapsed_serial.count() * 1000 << " milliseconds" << endl;
}

// Serial Merge Sort Execution
void serial_merge_sort_execution(const vector<int>& arr) {
    vector<int> serial_arr(arr);
    auto start_serial = chrono::steady_clock::now();
    merge_sort_serial(serial_arr, 0, serial_arr.size() - 1);
    auto end_serial = chrono::steady_clock::now();
    chrono::duration<double> elapsed_serial = end_serial - start_serial;
    // cout << "Serial Merge Sort Execution Time: " << elapsed_serial.count() * 1000 << " milliseconds" << endl;
}

// Serial Heap Sort Execution
void serial_heap_sort_execution(const vector<int>& arr) {
    vector<int> serial_arr(arr);
    auto start_serial = chrono::steady_clock::now();
    heap_sort_serial(serial_arr);
    auto end_serial = chrono::steady_clock::now();
    chrono::duration<double> elapsed_serial = end_serial - start_serial;
    // cout << "Serial Heap Sort Execution Time: " << elapsed_serial.count() * 1000 << " milliseconds" << endl;
}

// Parallel Bubble Sort Execution
void parallel_bubble_sort_execution(const vector<int>& arr) {
    int size = arr.size();
    int* d_arr;
    cudaMalloc(&d_arr, size * sizeof(int));
    cudaMemcpy(d_arr, arr.data(), size * sizeof(int), cudaMemcpyHostToDevice);
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    bubble_sort_kernel<<<1, size>>>(d_arr, size);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    // cout << "Parallel Bubble Sort Execution Time: " << milliseconds << " milliseconds" << endl;
    cudaFree(d_arr);
}

// Parallel Selection Sort Execution
void parallel_selection_sort_execution(const vector<int>& arr) {
    int size = arr.size();
    int* d_arr;
    cudaMalloc(&d_arr, size * sizeof(int));
    cudaMemcpy(d_arr, arr.data(), size * sizeof(int), cudaMemcpyHostToDevice);
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    selection_sort_kernel<<<1, size>>>(d_arr, size);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    // cout << "Parallel Selection Sort Execution Time: " << milliseconds << " milliseconds" << endl;
    cudaFree(d_arr);
}

// Parallel Insertion Sort Execution
void parallel_insertion_sort_execution(const vector<int>& arr) {
    int size = arr.size();
    int* d_arr;
    cudaMalloc(&d_arr, size * sizeof(int));
    cudaMemcpy(d_arr, arr.data(), size * sizeof(int), cudaMemcpyHostToDevice);
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    insertion_sort_kernel<<<1, size>>>(d_arr, size);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    // cout << "Parallel Insertion Sort Execution Time: " << milliseconds << " milliseconds" << endl;
    cudaFree(d_arr);
}

// Parallel Quick Sort Execution
void parallel_quick_sort_execution(const vector<int>& arr) {
    int size = arr.size();
    int* d_arr;
    cudaMalloc(&d_arr, size * sizeof(int));
    cudaMemcpy(d_arr, arr.data(), size * sizeof(int), cudaMemcpyHostToDevice);
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    parallel_quick_sort_kernel<<<1, size>>>(d_arr, size);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    // cout << "Parallel Quick Sort Execution Time: " << milliseconds << " milliseconds" << endl;
    cudaFree(d_arr);
}

// Parallel Merge Sort Execution
void parallel_merge_sort_execution(const vector<int>& arr) {
    int size = arr.size();
    int* d_arr;
    cudaMalloc(&d_arr, size * sizeof(int));
    cudaMemcpy(d_arr, arr.data(), size * sizeof(int), cudaMemcpyHostToDevice);
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    parallel_merge_sort_kernel<<<1, size>>>(d_arr, size);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    // cout << "Parallel Merge Sort Execution Time: " << milliseconds << " milliseconds" << endl;
    cudaFree(d_arr);
}

// Parallel Heap Sort Execution
void parallel_heap_sort_execution(const vector<int>& arr) {
    int size = arr.size();
    int* d_arr;
    cudaMalloc(&d_arr, size * sizeof(int));
    cudaMemcpy(d_arr, arr.data(), size * sizeof(int), cudaMemcpyHostToDevice);
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    parallel_heap_sort_kernel<<<1, size>>>(d_arr, size);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    // cout << "Parallel Heap Sort Execution Time: " << milliseconds << " milliseconds" << endl;
    cudaFree(d_arr);
}

int main() {
    int size;
    cout << "Enter the size of the array: ";
    cin >> size;

    vector<int> arr = generate_random_array(size);

    cout << "Algorithm\tSerial Time (ms)\tParallel Time (ms)\tSpeedup\t\tEfficiency" << endl;
    cout << "--------------------------------------------------------------------------------------------" << endl;

    // Bubble Sort
    auto start_serial_bubble = chrono::steady_clock::now();
    bubble_sort_serial(arr);
    auto end_serial_bubble = chrono::steady_clock::now();
    chrono::duration<double, milli> elapsed_serial_bubble = end_serial_bubble - start_serial_bubble;

    auto start_parallel_bubble = chrono::steady_clock::now();
    parallel_bubble_sort_execution(arr);
    auto end_parallel_bubble = chrono::steady_clock::now();
    chrono::duration<double, milli> elapsed_parallel_bubble = end_parallel_bubble - start_parallel_bubble;

    double speedup_bubble = elapsed_serial_bubble.count() / elapsed_parallel_bubble.count();
    double efficiency_bubble = speedup_bubble / 8; // Assuming 8 cores

    cout << "Bubble Sort\t" << fixed << setprecision(6) << setw(15) << elapsed_serial_bubble.count() << "\t   " << setw(15) << elapsed_parallel_bubble.count() << "\t\t" << setw(10) << speedup_bubble << "\t" << setw(10) << efficiency_bubble << endl;

    // Selection Sort
    auto start_serial_selection = chrono::steady_clock::now();
    selection_sort_serial(arr);
    auto end_serial_selection = chrono::steady_clock::now();
    chrono::duration<double, milli> elapsed_serial_selection = end_serial_selection - start_serial_selection;

    auto start_parallel_selection = chrono::steady_clock::now();
    parallel_selection_sort_execution(arr);
    auto end_parallel_selection = chrono::steady_clock::now();
    chrono::duration<double, milli> elapsed_parallel_selection = end_parallel_selection - start_parallel_selection;

    double speedup_selection = elapsed_serial_selection.count() / elapsed_parallel_selection.count();
    double efficiency_selection = speedup_selection / 8; // Assuming 8 cores

    cout << "Selection Sort\t" << fixed << setprecision(6) << setw(15) << elapsed_serial_selection.count() << "\t   " << setw(15) << elapsed_parallel_selection.count() << "\t\t" << setw(10) << speedup_selection << "\t" << setw(10) << efficiency_selection << endl;

    // Insertion Sort
    auto start_serial_insertion = chrono::steady_clock::now();
    insertion_sort_serial(arr);
    auto end_serial_insertion = chrono::steady_clock::now();
    chrono::duration<double, milli> elapsed_serial_insertion = end_serial_insertion - start_serial_insertion;

    auto start_parallel_insertion = chrono::steady_clock::now();
    parallel_insertion_sort_execution(arr);
    auto end_parallel_insertion = chrono::steady_clock::now();
    chrono::duration<double, milli> elapsed_parallel_insertion = end_parallel_insertion - start_parallel_insertion;

    double speedup_insertion = elapsed_serial_insertion.count() / elapsed_parallel_insertion.count();
    double efficiency_insertion = speedup_insertion / 8; // Assuming 8 cores

    cout << "Insertion Sort\t" << fixed << setprecision(6) << setw(15) << elapsed_serial_insertion.count() << "\t   " << setw(15) << elapsed_parallel_insertion.count() << "\t\t" << setw(10) << speedup_insertion << "\t" << setw(10) << efficiency_insertion << endl;

    // Quick Sort
    auto start_serial_quick = chrono::steady_clock::now();
    quick_sort_serial(arr, 0, arr.size() - 1);
    auto end_serial_quick = chrono::steady_clock::now();
    chrono::duration<double, milli> elapsed_serial_quick = end_serial_quick - start_serial_quick;

    auto start_parallel_quick = chrono::steady_clock::now();
    parallel_quick_sort_execution(arr);
    auto end_parallel_quick = chrono::steady_clock::now();
    chrono::duration<double, milli> elapsed_parallel_quick = end_parallel_quick - start_parallel_quick;

    double speedup_quick = elapsed_serial_quick.count() / elapsed_parallel_quick.count();
    double efficiency_quick = speedup_quick / 8; // Assuming 8 cores

    cout << "Quick Sort\t" << fixed << setprecision(6) << setw(15) << elapsed_serial_quick.count() << "\t   " << setw(15) << elapsed_parallel_quick.count() << "\t\t" << setw(10) << speedup_quick << "\t" << setw(10) << efficiency_quick << endl;

    // Merge Sort
    auto start_serial_merge = chrono::steady_clock::now();
    merge_sort_serial(arr, 0, arr.size() - 1);
    auto end_serial_merge = chrono::steady_clock::now();
    chrono::duration<double, milli> elapsed_serial_merge = end_serial_merge - start_serial_merge;

    auto start_parallel_merge = chrono::steady_clock::now();
    parallel_merge_sort_execution(arr);
    auto end_parallel_merge = chrono::steady_clock::now();
    chrono::duration<double, milli> elapsed_parallel_merge = end_parallel_merge - start_parallel_merge;

    double speedup_merge = elapsed_serial_merge.count() / elapsed_parallel_merge.count();
    double efficiency_merge = speedup_merge / 8; // Assuming 8 cores

    cout << "Merge Sort\t" << fixed << setprecision(6) << setw(15) << elapsed_serial_merge.count() << "\t   " << setw(15) << elapsed_parallel_merge.count() << "\t\t" << setw(10) << speedup_merge << "\t" << setw(10) << efficiency_merge << endl;

    // Heap Sort
    auto start_serial_heap = chrono::steady_clock::now();
    heap_sort_serial(arr);
    auto end_serial_heap = chrono::steady_clock::now();
    chrono::duration<double, milli> elapsed_serial_heap = end_serial_heap - start_serial_heap;

    auto start_parallel_heap = chrono::steady_clock::now();
    parallel_heap_sort_execution(arr);
    auto end_parallel_heap = chrono::steady_clock::now();
    chrono::duration<double, milli> elapsed_parallel_heap = end_parallel_heap - start_parallel_heap;

    double speedup_heap = elapsed_serial_heap.count() / elapsed_parallel_heap.count();
    double efficiency_heap = speedup_heap / 8; // Assuming 8 cores

    cout << "Heap Sort\t" << fixed << setprecision(6) << setw(15) << elapsed_serial_heap.count() << "\t   " << setw(15) << elapsed_parallel_heap.count() << "\t\t" << setw(10) << speedup_heap << "\t" << setw(10) << efficiency_heap << endl;

    return 0;
}
