/********************************************
 * 用于求解Lveine <<Quantum Chemistry>> 第7版
 * 14.3节中例题的程序，其作用为给定两个确定的
 * Slater型轨道，以此对氦原子做Hartree-Fock计
 * 算。
 * 计算结果与给出的结果完全一致。
 * 写于2021年10月13日。
********************************************/

# include <stdio.h>
# include <math.h>
# include <time.h>

double repulsion_int(int r, int s, int t, int u);
double overlap_int(int r, int s);
double Hcore_int(int r, int s);
double calc_fock(double F[][2], double P[][2]);
double solve_energy(double F[][2], double S[][2], double energy[2]);
double HF_energy(double P[][2], double *energy);

const double zeta1 = 1.45363;
const double zeta2 = 2.91093;

int main(void)
{

    double coeff[2] = {1, 1};
    double coeff_old[2] = {0, 0};
    double P[2][2];
    double F[2][2];
    double S[2][2];
    double energy[2]; 
    double k;

    int i1, i2;
    int count = 0;

    clock_t time1 = clock();

    //计算重叠矩阵
    for(i1=0; i1<2; i1++)
        for(i2=0; i2<2; i2++)
            S[i1][i2] = overlap_int(i1, i2);

    //给出一组初始系数，并归一化
    coeff[0] = 2, coeff[1] = 1;
    double norm = sqrt(coeff[0]*coeff[0] + coeff[1]*coeff[1] + 
                       2*coeff[0]*coeff[1]*overlap_int(1, 2));
    coeff[0] = coeff[0] / norm, coeff[1] = coeff[1] / norm;

    while(fabs(coeff[0] - coeff_old[0]) > 1e-8 && fabs(coeff[1] - coeff_old[1]) > 1e-8)
    {
        //备份原系数
        coeff_old[0] = coeff[0];
        coeff_old[1] = coeff[1];

        //计算P矩阵
        P[0][0] = 2 * coeff[0] * coeff[0];
        P[0][1] = 2 * coeff[0] * coeff[1];
        P[1][0] = 2 * coeff[1] * coeff[0];
        P[1][1] = 2 * coeff[1] * coeff[1];
    
        //计算Fock矩阵
        calc_fock(F, P);
    
        //求解轨道能量
        solve_energy(F, S, energy);
    
        //更新轨道系数
        k = - (F[1][0] - S[1][0] * energy[0]) / (F[1][1] - S[1][1] * energy[0]);
        //printf("%10lf\n", k);
        coeff[0] = 1 / sqrt(1 + k*k + 2*k*S[0][1]);
        coeff[1] = k / sqrt(1 + k*k + 2*k*S[0][1]);

        count ++;

        printf("%16.10lf%16.10lf\n", coeff[0], coeff[1]);
        printf("%16.10lf%16.10lf\n", coeff_old[0], coeff_old[1]);
        printf("%16.10lf%16.10lf\n", fabs(coeff[0] - coeff_old[0]), fabs(coeff[1] - coeff_old[1]));
        printf("E(HF) = %16.12lf\n\n", HF_energy(P, energy));
    }
    clock_t time2 = clock();
    printf("\nFinal coefficient: %16.10lf%16.10lf\n", coeff[0], coeff[1]);
    printf("cycles: %d\n", count);
    printf("Time usage: %.8lf ms\n", (double) (time2 - time1) / 1000);

    return 0;
}
//根据给出的公式计算排斥积分
double repulsion_int(int r, int s, int t, int u)
{
    const int judge = 1000 * r + 100 * s + 10 * t + u;
    double med1, med2, med3, med4, med5, med6;
    int i, j;
    if(judge == 1111)
        return 5 * zeta1 / 8;
    if(judge == 2222)
        return 5 * zeta2 / 8;
    if(judge == 1122 || judge == 2211)
    {
        med1 = zeta1 * zeta1 * zeta1 * zeta1 * zeta2;
        med2 = zeta1 * zeta1 * zeta1 * zeta2 * zeta2;
        med3 = zeta1 * zeta1 * zeta2 * zeta2 * zeta2;
        med4 = zeta1 * zeta2 * zeta2 * zeta2 * zeta2;
        med5 = 1;
        for(i=0; i<4; i++)
            med5 *= zeta1 + zeta2;
        return (med1 + med2 * 4 + med3 * 4 + med4) / med5;
    }
    if(judge == 1212 || judge == 2112 || judge == 1221 || judge == 2121)
    {
        med1 = zeta1 * zeta1 * zeta1;
        med2 = zeta2 * zeta2 * zeta2;
        med3 = 1;
        for(i = 0; i < 5; i++)
            med3 *= zeta1 + zeta2;
        return 20 * med1 * med2 / med3;
    }
    if(judge == 1112 || judge == 1121 || judge == 1211 || judge == 2111)
    {
        med1 = 16 * zeta1 * zeta1 * zeta1 * zeta1 * zeta2 * sqrt(zeta1) * sqrt(zeta2);
        med2 = 1;
        for(i=0; i<4; i++)
            med2 *= zeta1 * 3 + zeta2;
        med3 = 12 * zeta1 + 8 * zeta2;
        med4 = (zeta1 + zeta2) * (zeta1 + zeta2);
        med5 = 9 * zeta1 + zeta2;
        med6 = 2 * zeta1 * zeta1;
        return med1 / med2 *(med3 / med4 + med5 / med6);
    }
    if(judge == 1222 || judge == 2212 || judge == 2122 || judge == 2221)
    {
        med1 = 16 * zeta1 * zeta2 * zeta2 * zeta2 * zeta2 * sqrt(zeta1) * sqrt(zeta2);
        med2 = 1;
        for(i=0; i<4; i++)
            med2 *= zeta1 + zeta2 * 3;
        med3 = 8 * zeta1 + 12 * zeta2;
        med4 = (zeta1 + zeta2) * (zeta1 + zeta2);
        med5 = zeta1 + 9 * zeta2;
        med6 = 2 * zeta2 * zeta2;
        return med1 / med2 * (med3 / med4 + med5 / med6);
    }
}
//计算重叠积分
double overlap_int(int r, int s)
{
    if(r == s)
        return 1;
    else
    {
        double med1 = zeta1 * sqrt(zeta1) * zeta2 * sqrt(zeta2);
        double med2 = (zeta1 + zeta2) * (zeta1 + zeta2) * (zeta1 + zeta2);
        return 8 * med1 / med2;
    }
}
//计算H_{rs}^{core}
double Hcore_int(int r, int s)
{
    if(r == 1 && s == 1)
        return zeta1 * zeta1 / 2 - 2 * zeta1;
    else if(r == 2 && s == 2)
        return zeta2 * zeta2 / 2 - 2 * zeta2;
    else
    {
        double med1 = zeta1 * sqrt(zeta1) * zeta2 * sqrt(zeta2);
        double med2 = 4 * zeta1 * zeta2 - 8 * zeta1 - 8 * zeta2;
        double med3 = (zeta1 + zeta2) * (zeta1 + zeta2) * (zeta1 + zeta2);
        return med1 * med2 / med3;
    }
}
//计算Fock矩阵
double calc_fock(double F[][2], double P[][2])
{
    F[0][0] = Hcore_int(1, 1) + P[0][0] * repulsion_int(1, 1, 1, 1) / 2
              + P[0][1] * repulsion_int(1, 1, 1, 2) 
              + P[1][1] * (repulsion_int(1, 1, 2, 2) - repulsion_int(1, 2, 2, 1) / 2);

    F[0][1] = Hcore_int(1, 2) + P[0][0] * repulsion_int(1, 2, 1, 1) / 2
              + P[0][1] * (repulsion_int(1, 2, 2, 1) * 1.5 - repulsion_int(1, 1, 2, 2) / 2)
              + P[1][1] * repulsion_int(1, 2, 2, 2) / 2;

    F[1][0] = F[0][1];

    F[1][1] = Hcore_int(2, 2)
              + P[0][0] * (repulsion_int(2, 2, 1, 1) - repulsion_int(2, 1, 1, 2) / 2) 
              + P[0][1] * repulsion_int(2, 2, 1, 2) 
              + P[1][1] * repulsion_int(2, 2, 2, 2) / 2;
}
//求解轨道能量
double solve_energy(double F[][2], double S[][2], double energy[2])
{
    double a = 1 - S[0][1] * S[1][0];
    double b = F[1][0] * S[0][1] + F[0][1] * S[1][0] - F[0][0] - F[1][1];
    double c = F[0][0] * F[1][1] - F[0][1] * F[1][0];
    if(a > 0)
    {
        energy[0] = (-b - sqrt(b*b - 4 * a * c)) / (2 * a);
        energy[1] = (-b + sqrt(b*b - 4 * a * c)) / (2 * a);
    }
    else
    {
        energy[0] = (-b + sqrt(b*b - 4 * a * c)) / (2 * a);
        energy[1] = (-b - sqrt(b*b - 4 * a * c)) / (2 * a);
    }
}
//计算HF能量
double HF_energy(double P[][2], double *energy)
{
    double E_HF = 0;
    int i, j;
    for(i=0; i<2; i++)
        for(j=0; j<2; j++)
            E_HF += P[i][j] * Hcore_int(i+1, j+1) / 2;
    return E_HF + energy[0];
}
