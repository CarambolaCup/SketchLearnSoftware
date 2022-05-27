#include <iostream>
#include <algorithm>
#include <string.h>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <math.h>
#include <set>
#include <vector>

using namespace std;
#pragma warning(disable : 4996)
// 兰佳晨
// Encoded in CRLF UTF-8

// 流键 13个byte
// 时间戳 13个byte
const int ID_length = 13;
const int TimeStamp_length = 13;

class ID_input
{
public:
    char x[ID_length + TimeStamp_length];
    operator char* () const
    {
        return (char*)x;
    }
    operator unsigned char* () const
    {
        return (unsigned char*)x;
    }
};

// 未知何用,按字典序比较大小，和学姐提供的opeator相反
bool operator<(ID_input An, ID_input Bn)
{
    for (int i = 0; i < 26; i++)
    {
        if (An.x[i] < Bn.x[i])
        {
            return true;
        }
        else if (An.x[i] > Bn.x[i])
        {
            return false;
        }
    }
    return false;
}

// 将l,r,c参数及hash函数提前,以方便使用
const int l = 105;
const int r = 2;   // 2是我瞎写的
const int c = 100; // 100是我瞎写的
// r个哈希函数，需要搭框架的时候顺便实现一下(这里没有实现)
uint32_t test_hash_0(char* f)
{
    return (unsigned int)f[0] % (unsigned int)37;
}
uint32_t test_hash_1(char* f)
{
    return (unsigned int)f[0] % (unsigned int)57;
}
uint32_t (* hash_function[r + 1])(char*)={0,test_hash_0 ,test_hash_1 };

// 只是用来测试的两个哈希函数，37和57是随便选的

vector<ID_input> all_id_flow;
// 按V[k][i][j]排列 应该可以不用加1(但我还是加了)
unsigned int V[l + 1][r + 1][c + 1];
// 均值
double p[l + 1];
// 标准差,注意是方差开方
double sigma[l + 1];

//.dat 转 vector<ID>
int Read_Flowdata()
{
    char datafileName[100];
    // 注意文件路径
    sprintf(datafileName, "./0.dat");
    ID_input tmp_five_tuple;

    FILE* fin = fopen(datafileName, "rb");
    if (NULL != fin)
    {
        int k_count = 0;
        while (fread(&tmp_five_tuple, ID_length, 1, fin)) // 读13byte
        {
            all_id_flow.push_back(tmp_five_tuple);
            // 跳过时间戳
            fread(&tmp_five_tuple, TimeStamp_length, 1, fin);
            k_count++;
        }

        fclose(fin);
        // 约 2000 0000
        int test_a = all_id_flow.size();
        return 0;
    }
    else
    {
        return -1;
    }
}

int get_bit(unsigned char* a, int pos)
{
    int byte = pos / 8;
    int bit = pos % 8;
    if ((a[byte] & (1 << (7 - bit))) == 0)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

void Flow2Sketch()
{
    vector<ID_input> tmp_data = all_id_flow;
    for (ID_input tmp_flow : tmp_data)
    {
        for (size_t k = 0; k < ID_length * 8; k++)
        {
            // 0 则对 V[k] 无影响
            if (0 != get_bit((unsigned char*)tmp_flow, k))
            {
                for (size_t i = 0; i < r; i++)
                {
                    ++V[k][i][hash_function[i](tmp_flow)];
                }
            }
        }
    }
}

// 方差计算有问题，在此数据集下，由于数据过大会产生NaN 故加double损失精度
void Sketch2N_p_sigma()
{
    // 此处是否需要unsigned long long ?
    unsigned int sum;
    //平方和-和平方 未知精度是否足够
    double square_sum;
    for (size_t k = 0; k < ID_length * 8; k++)
    {
        sum = 0;
        square_sum = 0;
        for (size_t i = 0; i < r; i++)
        {
            for (size_t j = 0; j < c; j++)
            {
                sum += V[k][i][j];
                square_sum += (double)V[k][i][j] * (double)V[k][i][j]; // 约10^12数量级
            }
        }
        p[k] = (double)sum / (double)(r * c);
        sigma[k] = sqrt(square_sum / (long double)(r * c) - p[k] * p[k]);
    }
}

// 杨仕博：代码里bit_flow、prob_vector、每bit概率是从1开始的，以契合算法中的表述
//         flow是从0开始的，以契合给定的数据和哈希函数

double normalCFD(double value) // 计算标准正态分布
{
    return 0.5 * erfc(-value / sqrt(2));
}

struct ans_t // 是extract large flow返回向量中元素的type
{
    char bit_flow[l + 2];       // 这个是用1个byte的'0'或'1'表示1个bit
    char flow[(l + 7) / 8 + 1]; // 这个是跟源数据相同的结构
    unsigned int size;          // 流的大小
    double prob_vector[l + 2];  // 可能性向量
    ans_t(char* bbit_flow, char* fflow, unsigned int ssize = 0, double* pprob_vector = NULL)
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
    unsigned int V[][r+1][c+1], double* p, double* sigama2, int k) // 计算大流第k个bit为1的概率
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

struct two_types_of_flow // 两种流表示
{
    char bit_flow[l + 2];       // 一个'1'或者'0'表示一个bit
    char flow[(l + 7) / 8 + 1]; // 应该是和给定的数据一样的结构(不知道会不会有bug...)
    two_types_of_flow(char* bbit_flow, char* fflow)
    {
        strcpy(bbit_flow, bit_flow);
        strcpy(fflow, flow);
    }
};

vector<two_types_of_flow> possible_flows; // 一个辅助函数的全局变量
char current_T[l + 2];                    // 一个辅助函数的全局变量

int flag_ = 0;

void find_possible_flows(int i, int j, int k, char* T) // 找到正则表达式中所有可能的流
{
    if (k == l + 1)
    {
        char ans[(l + 7) / 8 + 1];
        
        for (int ii = 0; ii < (l + 7) / 8; ii++)
        {
            ans[ii] = 0;
            for (int jj = 0; jj < 8; jj++)
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
            if (i == 1 && j == 57)printf("ii=%d completed\n", ii);
        }
        ans[(l + 7) / 8] = '\0';
        if (i == 1 && j == 57)printf("V[%d][%d] before hash\n", i, j);
        if (hash_function[i](ans) % c == j)
        {
            possible_flows.push_back(two_types_of_flow(current_T, ans));
        }
        return;
    }
    else
    {
        if (T[k] != '*')
        {
            current_T[k] = T[k];
            if (i == 1 && j == 57)printf("V[%d][%d] %d completed,is%c\n", i, j,k,current_T[k]);
            find_possible_flows(i, j, k + 1, T);
        }
        else
        {
            current_T[k] = '0';
            if (i == 1 && j == 57)printf("V[%d][%d] %d is *, = %c\n", i, j, k, current_T[k]);
            find_possible_flows(i, j, k + 1, T);
            current_T[k] = '1';
            if (i == 1 && j == 57)printf("V[%d][%d] %d is *, = %c\n", i, j, k, current_T[k]);
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
    unsigned int V[][r+1][c+1], double* p, double* sigama2)
{
    // 第一步，计算每个bit的概率估值
    double hat_p[l+1];
    for (int k = 1; k <= l; k++)
    {
        hat_p[k] = cal_hat_p(theta, i, j, V, p, sigama2, k);
    }
    printf("V[%d][%d] step 1 completed\n", i, j);
    // 第二步，找到所有候选的大流，存在possible_flows里面
    char T[l + 2];
    for (int k = 1; k <= l; k++)
    {
        if (hat_p[k] > 0.99)
        {
            T[k] = '1';
        }
        else if (1 - hat_p[k] > 0.99)
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
    if(i==1&&j==57)printf("V[%d][%d] before find_possible_flows\n", i, j);
    find_possible_flows(i, j, 1, T);
    printf("V[%d][%d] step 2 completed\n", i, j);
    // 第三步，估计大流的频率和可能性向量
    double estimated_frequency[l + 1];
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
                estimated_frequency[k] = (1 - r / p[k]) * V[0][i][j];
                estimated_p[k] = 1 - hat_p[k];
            }
        }
        sort(estimated_frequency + 1, estimated_frequency + 1 + l);
        result.push_back(ans_t(item->bit_flow, item->flow, estimated_frequency[l / 2], estimated_p));
    }
    printf("V[%d][%d] step 3 completed\n", i, j);
    // 第四步，去sketch里查候选流的数据，删掉过小的
    if (r == 1)
    {
        return result;
    }
    for (vector<ans_t>::iterator item = result.begin(); item != result.end(); item++)
    {
        if (flag_ == 1)
        {
            vector<ans_t>::iterator tmp_item = item;
            tmp_item--;
            result.erase(tmp_item);
        }
        flag_ = 0;
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
            flag_ = 1;
            /*
            vector<ans_t>::iterator tmp_item = item;
            item--;
            result.erase(tmp_item);
            */
        }
    }
    return result;
}

//从sketch中删去大流
void RemoveFlows(vector<ans_t> FF) {
    int FF_size = FF.size();
    for (int it = 0; it < FF.size(); it++) {
        char* ans = FF[it].flow;
        for (size_t k = 0; k < ID_length * 8; k++)
        {
            if (0 != get_bit((unsigned char*)ans, k))
            {
                for (size_t i = 0; i < r; i++)
                {
                    int j = hash_function[i](ans);
                    V[k][i][j] -= FF[it].size;
                }
            }
        }
    } 
}
//检查当前sketch是否符合高斯分布
//我理解的是所有L个sketch每一个都符合1sigma,2sigma,3sigma的数量要求
bool Terminate() {
    for (size_t k = 0; k < ID_length * 8; k++) {
        size_t sigma_num1 = 0, sigma_num2 = 0, sigma_num3 = 0;
        for (int i = 1; i <= r; i++) {
            for (int j = 1; j <= c; j++) {
                if (V[k][i][j] <= p[k] + 3.0 * sigma[k] && V[k][i][j] >= p[k] - 3.0 * sigma[k])sigma_num3++;
                if (V[k][i][j] <= p[k] + 2.0 * sigma[k] && V[k][i][j] >= p[k] - 2.0 * sigma[k])sigma_num2++;
                if (V[k][i][j] <= p[k] + 1.0 * sigma[k] && V[k][i][j] >= p[k] - 1.0 * sigma[k])sigma_num1++;
            }
        }
        if ((double)sigma_num1 / r * c < 0.6826)return false;
        if ((double)sigma_num2 / r * c < 0.9544)return false;
        if ((double)sigma_num3 / r * c < 0.9973)return false;
    }
    return true;
}
int main()
{
    if (0 == Read_Flowdata()) //流数据读入
    {
        // 赋予哈希函数
        hash_function[0] = test_hash_0;
        hash_function[1] = test_hash_1;
        // 流转sketch
        Flow2Sketch();
        // sketch 生成N(p,sigma)
        Sketch2N_p_sigma();
        printf("Read in flow data success\n");
    }
    else
    {
        printf("Read in flow data error\n");
    }

    //F为所有大流集合，FF为每次循环找出的大流集合
    vector<ans_t> F;
    double theta = 0.5;
    int nnnn = 4;
    //while (nnnn--) {
        vector<ans_t> FF;
        for (int i = 1; i <= r; i++) {
            for (int j = 1; j <= c; j++) {
                vector<ans_t> temp_F = ExtractLargeFlows(theta, i, j,
                    V, p, sigma);
                if (!temp_F.empty()) {
                    for (vector<ans_t>::iterator it = temp_F.begin(); it < temp_F.end(); it++)
                        FF.push_back(*it);
                }
                printf("V[%d][%d]completed\n",i,j);
            }
        }
        //本次循环找出大流时，剔除大流，重新计算期望、方差
        if (!FF.empty()) {
            for (vector<ans_t>::iterator it = FF.begin(); it < FF.end(); it++)
                F.push_back(*it);
            RemoveFlows(FF);
            printf("RemoveFlowscompleted\n");
            Sketch2N_p_sigma();
        }
       // if (Terminate())break;
        //没有找出大流，theta减半
        //if (FF.empty())theta /= 2;
    //}
    return 0;
}