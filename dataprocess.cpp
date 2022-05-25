#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <set>

using namespace std;

class ID_input
{
public:
    char x[26] = {0};
};

bool operator<(ID_input an, ID_input bn)
{
    for (int i = 0; i < 26; i++)
    {
        if (bn.x[i] < an.x[i])
        {
            return true;
        }
        else if (bn.x[i] > an.x[i])
        {
            return false;
        }
    }
    return false;
}

int main()
{
    set<ID_input> all_id_set;
    char datafileName[100];
    sprintf(datafileName, "./data/0.dat");
    FILE *fin = fopen(datafileName, "rb");
    ID_input tmp_five_tuple;

    int kkk = 0;
    while (fread(&tmp_five_tuple, 1, 8, fin))
    {
        if (kkk % 2 == 1)
        {
            all_id_set.insert(tmp_five_tuple);
        }
        ++kkk;
    }
    int test_a = all_id_set.size();
        
    fclose(fin);
    return 0;
}