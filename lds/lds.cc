#include<bits/stdc++.h>
using namespace std;

int LDS(int a[],int size){
    int dp[size+1],res=1;
    for(int i=0;i<size;i++){
        dp[i] = 1;
    }
    for(int i=size-1;i>=0;i--){
        for(int j=0;j<i;j++){
            if(a[i]>a[j]){
                dp[i] = max(dp[i],1+dp[j]);
            }
        }
        res = max(res,dp[i]);
    }

    for(int i=0;i<size;i++){
        cout<<dp[i]<<" ";
    }
    cout<<endl;
    
    return res;
}

int LDS2(int a[],int size){
    int dp[size+1],res=1;
    for(int i=0;i<size;i++){
        dp[i] = 1;
    }
    for (int c0 = 1; c0 < size; c0 += 1)
        for (int c1 = size - c0; c1 < size; c1 += 1)
        {
            if(a[c1]>a[size - c0 - 1]){
                dp[c1] = max(dp[c1], 1+dp[size - c0 - 1]);
            }
            //(c1, size - c0 - 1);
        }
    for(int i=0;i<size;i++){
        cout<<dp[i]<<" ";
    }
    cout<<endl;
    
    return res;
}

int main(){
    int a[] = {10,9,2,5,3,7,101,18, 11, 15, 13};
    int size = sizeof(a)/sizeof(a[0]);
    cout<<LDS(a,size)<<"\n";
    cout<<LDS2(a,size)<<"\n";
}
