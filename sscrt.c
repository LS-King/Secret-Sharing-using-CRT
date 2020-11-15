#include <stdio.h>
#include "miracl.h"
//(3,5)����
#define N 5
#define T 3

void mycrt(big* d, big* k, int* recIndex, int equationNum, big testsecret);

int main()
{
	FILE* fp;
	int s_num, d_num;
	int equationNum = 0;
	int recIndex[N] = {0}; //���ڻָ��Ĳ��ظ��������ܶԱ��
	char sstr[1000];
	miracl* mip = mirsys(500, 0);
	//������ʼ��
	big secret = mirvar(0);//ԭʼ����
	big testsecret = mirvar(0);//�ָ�����
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
	cinnum(secret, fp);	//���ļ��ж�������s
	printf("ԭʼ����Ϊ��\n");
	cotnum(secret, stdout);
	s_num = cotstr(secret, sstr);//s��λ��
	d_num = s_num * 2 / (2 * T - 1);//һ�������d��λ��
	//���������
	irand((unsigned)time(NULL));
	for (int i = 0; i < N;){
		bigdig(d_num, 10, temp);
		if (isprime(temp)){
			copy(temp, d[i]);
			i++;
		}
	}
	//����
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
	//���������ܶ�
	for (int i = 0; i < N; i++){
		powmod(secret, one, d[i], k[i]);//secret = ki (mod di)
	}
	for (int i = 0; i < N; i++){
		printf("\n��%d�������ܣ�\n", i + 1);
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
		printf("\n���������ڻָ��������ܶԸ�����\n");
		scanf("%d", &equationNum);
		if(equationNum > N || equationNum <= 0){
			printf("��ֵ�Ƿ���\n");
		}
		else{
			break;
		}
	}
	printf("\n���������ڻָ��Ĳ��ظ��������ܶԱ�ţ�\n");
	int i = 0;
	while (i < equationNum){
		scanf("%d", &recIndex[i]);
		if(recIndex[i] > N || recIndex[i] <= 0){
			printf("�Ƿ����룡\n");
			continue;
		}
		else{
			recIndex[i] = recIndex[i] - 1;
			i++;
		}
	}
    
	mycrt(d, k, recIndex, equationNum, testsecret);//�й�ʣ�ඨ���ͬ�෽����
	printf("\n�ָ������\n");
	cotnum(testsecret, stdout);
	if (mr_compare(testsecret, secret) == 0){
		printf("\n�ɹ��ָ�����\n");
	}
	else{
		printf("\n���ָܻ�ʧ��\n");
	}
	mirexit();
	return 0;
}
/* �й�ʣ�ඨ�� *
 * k = secret (mod d) *
 * recIndex - ���ڻָ��Ĳ��ظ��������ܶԱ�� *
 * equationNum - ͬ�෽�̸��� *
 * testsecret - ��Կ�ָ���� */
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