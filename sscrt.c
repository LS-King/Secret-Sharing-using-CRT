#include <stdio.h>
#include "miracl.h"
//(3,5)门限
#define N 5
#define T 3

void mycrt(big* d, big* k, int* recIndex, int equationNum, big testsecret);

int main()
{
	FILE* fp;
	int s_num, d_num;
	int equationNum = 0;
	int recIndex[N] = {0}; //用于恢复的不重复的子秘密对标号
	char sstr[1000];
	miracl* mip = mirsys(500, 0);
	//变量初始化
	big secret = mirvar(0);//原始秘密
	big testsecret = mirvar(0);//恢复秘密
	big temp = mirvar(0);
	big one = mirvar(1);
	//n > k > m
	big n = mirvar(1);  
	big m = mirvar(1);
	big k[N], d[N];
	for (int i = 0; i < N; i++){
		k[i] = mirvar(0);
		d[i] = mirvar(0);
	}

	fp = fopen("secret.txt", "r+");
	cinnum(secret, fp);	//从文件中读入秘密s
	printf("原始秘密为：\n");
	cotnum(secret, stdout);
	s_num = cotstr(secret, sstr);//s的位数
	d_num = s_num * 2 / (2 * T - 1);//一个合理的d的位数
	//生成随机数
	irand((unsigned)time(NULL));
	for (int i = 0; i < N;){
		bigdig(d_num, 10, temp);
		if (isprime(temp)){
			copy(temp, d[i]);
			i++;
		}
	}
	//排序
	for (int i = 0; i < N - 1; i++){
		for (int j = i + 1; j < N; j++){
			if (mr_compare(d[i], d[j]) > 0){
				copy(d[i], temp);
				copy(d[j], d[i]);
				copy(temp, d[j]);
			}
		}
	}
	//n = d1 * d2 * ... * dt
	for (int i = 0; i < T; i++){
		multiply(n, d[i], n);
	}
	//m = d(n-t+2) * d(n-t+3) * ... * dn
	for (int i = N - T + 1; i < N; i++){
		multiply(m, d[i], m);
	}
	//生成子秘密对
	for (int i = 0; i < N; i++){
		powmod(secret, one, d[i], k[i]);//secret = ki (mod di)
	}
	for (int i = 0; i < N; i++){
		printf("\n第%d组子秘密：\n", i + 1);
		printf("d%d = ", i + 1);
		cotnum(d[i], stdout);
		printf("k%d = ", i + 1);
		cotnum(k[i], stdout);
	}

	printf("\nN = ");
	cotnum(n, stdout);
	printf("\nM = ");
	cotnum(m, stdout);

	while(1){
		printf("\n请输入用于恢复的子秘密对个数：\n");
		scanf("%d", &equationNum);
		if(equationNum > N || equationNum <= 0){
			printf("数值非法！\n");
		}
		else{
			break;
		}
	}
	printf("\n请输入用于恢复的不重复的子秘密对标号：\n");
	int i = 0;
	while (i < equationNum){
		scanf("%d", &recIndex[i]);
		if(recIndex[i] > N || recIndex[i] <= 0){
			printf("非法输入！\n");
			continue;
		}
		else{
			recIndex[i] = recIndex[i] - 1;
			i++;
		}
	}
    
	mycrt(d, k, recIndex, equationNum, testsecret);//中国剩余定理解同余方程组
	printf("\n恢复结果：\n");
	cotnum(testsecret, stdout);
	if (mr_compare(testsecret, secret) == 0){
		printf("\n成功恢复秘密\n");
	}
	else{
		printf("\n秘密恢复失败\n");
	}
	mirexit();
	return 0;
}
/* 中国剩余定理 *
 * k = secret (mod d) *
 * recIndex - 用于恢复的不重复的子秘密对编号 *
 * equationNum - 同余方程个数 *
 * testsecret - 密钥恢复结果 */
void mycrt(big* d, big* k, int* recIndex, int equationNum, big testsecret){
	big M[N], MI[N], x[N];
	for (int i = 0; i < N; i++)
	{
		M[i] = mirvar(0);
		MI[i] = mirvar(0);
		x[i] = mirvar(1);
	}
	big mul = mirvar(1);
	big temp = mirvar(0);
	big sum = mirvar(0);
	big one = mirvar(1);

	for (int i = 0; i < equationNum; i++){
		multiply(mul, d[recIndex[i]], mul);
	}
	copy(mul, temp);
	for (int i = 0; i < equationNum; i++){
		divide(mul, d[recIndex[i]], M[i]);
		copy(temp, mul);
	}
	for (int i = 0; i < equationNum; i++){
		xgcd(M[i], d[recIndex[i]], MI[i], MI[i], MI[i]);
	}
	for (int i = 0; i < equationNum; i++){
		multiply(x[i], M[i], x[i]);
		multiply(x[i], MI[i], x[i]);
		multiply(x[i], k[recIndex[i]], x[i]);
	}
	for (int i = 0; i < equationNum; i++){
		add(sum, x[i], sum);
	}
	powmod(sum, one, mul, sum);
	copy(sum, testsecret);
}