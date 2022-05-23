#include <iostream>
#include <algorithm>
#include <string.h>
#include <vector>
#include <cmath>
using namespace std;

// 杨仕博：代码里bit_flow、prob_vector、每bit概率是从1开始的，以契合算法中的表述
//         flow是从0开始的，以契合给定的数据和哈希函数

const int l = 105;
const int r = 2;            // 2是我瞎写的
const int c = 100;          // 100是我瞎写的

uint32_t(*hash_function[r])(char*);//r个哈希函数，需要搭框架的时候顺便实现一下(这里没有实现)

double normalCFD(double value)// 计算标准正态分布
{
    return 0.5 * erfc(-value / sqrt(2));
}

struct ans_t// 是extract large flow返回向量中元素的type
{
    char bit_flow[l+2];                 // 这个是用1个byte的'0'或'1'表示1个bit
    char flow[(l + 7) / 8 + 1];         // 这个是跟源数据相同的结构
    unsigned int size;                  // 流的大小
    double prob_vector[l+2];            // 可能性向量
    ans_t(char* bbit_flow,char* fflow,unsigned int ssize = 0, double* pprob_vector = NULL)
    {
        strcpy(bbit_flow, bit_flow);
        strcpy(fflow, flow);
        if (pprob_vector != NULL)
        {
            for (int i = 1; i <= l; i++)
            {
                prob_vector[i] = pprob_vector[i];
            }
        }
        size = ssize;
    }
};

double cal_hat_p(double theta, int i, int j,
    int*** V, double* p, double* sigama2, int k)// 计算大流第k个bit为1的概率
{
    double r = (double)V[k][i][j] / V[0][i][j];
    if (r < theta)
    {
        return 0;
    }
    if (1 - r < theta)
    {
        return 1;
    }
    double ans = 0;
    double prob_1 = (V[k][i][j] - theta * V[0][i][j]) /
        (V[k][i][j] - theta * V[0][i][j]);
    double prob_0 = (V[k][i][j]) /
        (V[k][i][j] - theta * V[0][i][j]);
    double normal_val1 = normalCFD((prob_1 - p[k]) / sqrt(sigama2[k]));
    double normal_val0 = normalCFD((prob_0 - p[k]) / sqrt(sigama2[k]));
    return normal_val1 * p[k] + (1 - normal_val0) * (1 - p[k]);
}

struct two_types_of_flow// 两种流表示
{
    char bit_flow[l + 2];                   // 一个'1'或者'0'表示一个bit
    char flow[(l + 7) / 8 + 1];             // 应该是和给定的数据一样的结构(不知道会不会有bug...)
    two_types_of_flow(char* bbit_flow, char* fflow)
    {
        strcpy(bbit_flow, bit_flow);
        strcpy(fflow, flow);
    }
};

vector<two_types_of_flow> possible_flows;   // 一个辅助函数的全局变量
char current_T[l + 2];                      // 一个辅助函数的全局变量

void find_possible_flows(int i, int j, int k, char* T)// 找到正则表达式中所有可能的流
{
    if (k == l + 1)
    {
        char ans[(l+7)/8 + 1];
        for (int ii = 0; ii < (l + 7) / 8; ii++)
        {
            ans[ii] = 0;
            for (int jj = 0; jj < 8; j++)
            {
                if (current_T[8 * ii + jj + 1] == '1')
                {
                    ans[ii] |= (1 << jj);
                }
                else
                {
                    ans[ii] &= (~(1 << jj));
                }
            }
        }
        ans[(l + 7) / 8] = '\0';
        if (hash_function[i](ans) % c == j)
        {
            possible_flows.push_back(two_types_of_flow(current_T ,ans));
        }
        return;
    }
    else
    {
        if (T[k] != '*')
        {
            current_T[i] = T[i];
            find_possible_flows(i, j, k + 1, T);
        }
        else
        {
            current_T[i] = '0';
            find_possible_flows(i, j, k + 1, T);
            current_T[i] = '1';
            find_possible_flows(i, j, k + 1, T);
        }
    }
    return;
}

/* 返回ans_t的一个向量，每个元素存储了：
 *      
 *      用byte表示bit的流键
 *      和输入数据同结构的流键
 *      流的大小
 *      可能性向量
 * 
 */
vector<ans_t> ExtractLargeFlows(double theta, int i, int j,
    int*** V, double* p, double* sigama2)
{
    // 第一步，计算每个bit的概率估值
    double hat_p[l];
    for (int k = 1; k <= l; k++)
    {
        hat_p[k] = cal_hat_p(theta, i, j, V, p, sigama2, k);
    }

    // 第二步，找到所有候选的大流，存在possible_flows里面
    char T[l+2];
    for (int k = 1; k <= l; k++)
    {
        if (hat_p[k] > 0.99)
        {
            T[k] = '1';
        }
        else if (1 - hat_p[k] < 0.99)
        {
            T[k] = '0';
        }
        else
        {
            T[k] = '*';
        }
    }
    T[l + 1] = '\0';
    current_T[l + 1] = '\0';
    possible_flows.clear();
    find_possible_flows(i, j, 1, T);

    // 第三步，估计大流的频率和可能性向量
    double estimated_frequency[l+1];
    double estimated_p[l + 1];
    vector<ans_t> result;
    for (vector<two_types_of_flow>::iterator item = possible_flows.begin();
        item != possible_flows.end(); item++)
    {
        for (int k = 1; k <= l; k++)
        {
            if (item->bit_flow[k] == '1')
            {
                double r = (double)V[k][i][j] / V[0][i][j];
                estimated_frequency[k] = ((r - p[k]) / (1 - p[k])) * V[0][i][j];
                estimated_p[k] = hat_p[k];
            }
            else
            {
                double r = (double)V[k][i][j] / V[0][i][j];
                estimated_frequency[k] = (1 - r/p[k]) * V[0][i][j];
                estimated_p[k] = 1 - hat_p[k];
            }
        }
        sort(estimated_frequency + 1, estimated_frequency + 1 + l);
        result.push_back(ans_t(item->bit_flow, item->flow, estimated_frequency[l / 2], estimated_p));
    }

    // 第四步，去sketch里查候选流的数据，删掉过小的
    if (r == 1)
    {
        return result;
    }
    for (vector<ans_t>::iterator item = result.begin(); item != result.end(); item++)
    {
        for (int ii = 1; ii <= r; ii++)
        {
            if (ii == i)
            {
                continue;
            }
            int jj = hash_function[ii](item->flow) % c;
            for (int k = 1; k <= l; k++)
            {
                if (item->bit_flow[k] == '0' && V[0][ii][jj] - V[k][ii][jj] < item->size)
                {
                    item->size = V[0][ii][jj] - V[k][ii][jj];
                }
                else if (item->bit_flow[k] == '1' && V[k][ii][jj] < item->size)
                {
                    item->size = V[k][ii][jj];
                }
            }
        }
        if (item->size < theta * V[0][i][j])
        {
            vector<ans_t>::iterator tmp_item = item;
            item--;
            result.erase(tmp_item);
        }
    }
    return result;
}

int main()
{

}