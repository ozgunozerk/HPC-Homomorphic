#include <iostream>
#include <omp.h>

using namespace std;


void sequential_intt(int &* f, int & N, int & Ninv, int & q, int &* zetas_inv){
    int j = 0;
    int k = 0;
    int N_ = N/2;
    int start;
    for(int length = 2; length <= N_; length << 1){
        int start = 0;
        for(start = 0; start < N; start + (j+1) + length){
            int omega = zetas_inv[k++];
            for(j = start; j < start+length; j++){
                int t = f[j];
                f[j] = (t+f[j + length])%q;
                f[j + length] = (t - f[j + length])%q;
                f[j + length] = (omega*f[j + length])%q;
            }
        }
    }
    for (int a=0; a < N; a++){
        f[j] = (f[j]*Ninv)%q; 
    }    
} 



void sequential_ntt(int** f, const int N, const int q, int** zetas)
{
//    cout << N << endl;
    for (int length = N / 2; length >= 2; length /= 2)
    {
//        cout << length << " | ";

        int k = N / (2 * length);
        for (int start = 0; start < N; start += length * 2)
        {
//            cout << start << ","<< k << " ";
            int omega = zetas[k][0];
            k += 1;

            for (int j = start; j < start + length; j++)
            {
                int t = (omega * f[j + length][0]) % q;
                f[j + length][0] = (f[j][0] - t + q) % q;
                f[j][0] = (f[j][0] + t) % q;
            }
        }

//        cout << endl;
    }
}

void parallel_ntt(int** f, const int N, const int q, int** zetas)
{
    for (int length = N / 2; length >= 2; length /= 2)
    {
        const int k_start = N / (length * 2);

#pragma omp parallel for default(shared)
        for (int i = 0; i < k_start; i++)
        {
            const int omega = zetas[k_start + i][0];
            const int f_start = i * length * 2;
            const int f_end = f_start + length;

            for (int j = f_start; j < f_end; j++)
            {
                int t = (omega * f[j + length][0]) % q;
                f[j + length][0] = (f[j][0] - t + q) % q;
                f[j][0] = (f[j][0] + t) % q;
            }
        }
    }
}






int main()
{
    omp_set_num_threads(8);
    const int N = 1u << 25u;
    const int q = 7681;
    const int padding = 40;
    int ** f_s = new int*[N];
    int ** f_p = new int*[N];
    int ** zetas = new int*[N];
    for (int a = 0; a < N; ++a)
    {
        f_s[a] = new int[padding];
        f_p[a] = new int[padding];
        zetas[a] = new int[padding];
    }

    for (int i = 0; i < N; ++i)
        f_s[i][0] = f_p[i][0] = zetas[i][0] = i;

    double start;

    start = omp_get_wtime();
    sequential_ntt(f_s, N, q, zetas);
    auto sequential_time = omp_get_wtime() - start;
    cout << "Sequential: " << sequential_time << " seconds" << endl;

    start = omp_get_wtime();
    parallel_ntt(f_p, N, q, zetas);
    auto parallel_time = omp_get_wtime() - start;
    cout << "Parallel: " << parallel_time << " seconds" << endl;

    auto speedup = sequential_time / parallel_time;
    cout << "Speed-up: " << speedup << endl;
    auto thread_count = omp_get_max_threads();
    auto efficiency = speedup / thread_count;
    cout << "Efficiency for " << thread_count << " threads: " << efficiency << endl;

    cout << "f_p == f_s: ";
    bool is_same = true;
    for (int i = 0; i < N && is_same; i++)
        is_same = (f_p[i][0] == f_s[i][0]);
    cout << (is_same ? "True" : "FALSE") << endl;

    for (int i = 0; i < N; ++i){
        delete[] f_s[i];
        delete[] f_p[i];
        delete[] zetas[i];
    }

    delete[] f_s;
    delete[] f_p;
    delete[] zetas;

    return 0;
}
