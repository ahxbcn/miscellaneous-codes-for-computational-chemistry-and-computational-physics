#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#define LENGTH 75 
#define RELAXCYCLE (LENGTH * LENGTH * 100)
#define RUNCYCLE (LENGTH * LENGTH * 1000)
#define PROP_PER_FRAME 1000

const int prmax = RAND_MAX - RAND_MAX % LENGTH;
const double T_min = 0.01;
int x, y;
double delta_E;
double T = 4.00;
short field[LENGTH * LENGTH] = {0};

int GetRandomInt(int range);
double GetTotalEner(void);
double GetTotalMagMom(void);
double GetRandomDoub(void);
void runstep(void);

int main(void)
{
    int i, j;
    int x, y;
    int count[LENGTH] = {0};
    int balanced_flag = 0;
    double energy, C, chi;
    double average_E, average_E2;
    double mag_moment, mag_moment2;
    double ave_mag_moment, ave_mag_moment2;
    double per_mag_moment, per_mag_moment2;
    double per_average_E, per_average_E2;

    double old_per_average_E;
    double old_per_ave_mag_moment;

    time_t time1, time2;
    register int randnumber = RAND_MAX;
    FILE *fp, *fp_E, *fp_m, *fp_C, *fp_chi, *fp_conf;

    time1 = clock();
    for(i=0; i<LENGTH*LENGTH; i++)
    {
        field[i] = rand() % 2 * 2 - 1;
    }
    fp = fopen("ising_model.txt", "w");
    fp_E = fopen("ising_model_E.txt", "w");
    fp_C = fopen("ising_model_C.txt", "w");
    fp_m = fopen("ising_model_m.txt", "w");
    fp_chi = fopen("ising_model_chi.txt", "w");
    while(T >= T_min)
    {
        printf("T = %lf\n", T);
        fprintf(fp, "T = %lf\n", T);
        printf("Relaxing\n");

        average_E = 0;
        ave_mag_moment = 0;
        old_per_average_E = GetTotalEner() / (LENGTH * LENGTH);
        old_per_ave_mag_moment = GetTotalMagMom() / (LENGTH * LENGTH);
        balanced_flag = 0;
        i = 0;

        srand(clock() + time(NULL));

        for(i=0; i<RELAXCYCLE; i++)
        {
            runstep();
            if(i % (LENGTH * LENGTH) == 0)
                printf("%8d", i);
            i++;
        }

        average_E = 0;
        average_E2 = 0;
        ave_mag_moment = 0;
        ave_mag_moment2 = 0;
        printf("\nProducting\n");
        for(i=0; i<RUNCYCLE; i++)
        {
            runstep();
            if(i % PROP_PER_FRAME == 0)
            {
                energy = GetTotalEner();
                average_E += energy * PROP_PER_FRAME / RUNCYCLE;
                average_E2 += energy * energy * PROP_PER_FRAME / RUNCYCLE;
                mag_moment = GetTotalMagMom();
                ave_mag_moment += mag_moment * PROP_PER_FRAME / RUNCYCLE;
                ave_mag_moment2 += mag_moment * mag_moment * PROP_PER_FRAME / RUNCYCLE;
            }
            if(i % (LENGTH * LENGTH * 10) == 0)
                printf("%8d", i);
        }
        per_average_E = average_E / (LENGTH * LENGTH);
        per_average_E2 = average_E2 / (LENGTH * LENGTH);
        C = (average_E2 - average_E * average_E) / (T * T * LENGTH * LENGTH);
        per_mag_moment = ave_mag_moment / (LENGTH * LENGTH);
        per_mag_moment2 = ave_mag_moment2 / (LENGTH * LENGTH);
        chi = (ave_mag_moment2 - ave_mag_moment * ave_mag_moment) / (T * LENGTH * LENGTH);
        printf("\n<E>: %lf\n", average_E);
        printf("<E2>: %lf\n", average_E2);
        printf("<E>/N: %lf\n", per_average_E);
        printf("<E2>/N: %lf\n", per_average_E2);
        printf("C: %lf\n", C);
        printf("<m>/N: %lf\n", ave_mag_moment / (LENGTH * LENGTH));
        printf("<m2>/N: %lf\n", ave_mag_moment2 / (LENGTH * LENGTH));
        printf("chi: %lf\n", chi);
        printf("\n");
        fprintf(fp, "<E>: %lf\n", average_E);
        fprintf(fp, "<E2>: %lf\n", average_E2);
        fprintf(fp, "<E>/N: %lf\n", per_average_E);
        fprintf(fp, "<E2>/N: %lf\n", per_average_E2);
        fprintf(fp, "C: %lf\n", C);
        fprintf(fp, "<m>/N: %lf\n", ave_mag_moment / (LENGTH * LENGTH));
        fprintf(fp, "<m2>/N: %lf\n", ave_mag_moment2 / (LENGTH * LENGTH));
        fprintf(fp, "chi: %lf\n", chi);
        fprintf(fp, "\n");

        fprintf(fp_E, "%6.4lf%20.12f\n", T, per_average_E);
        fprintf(fp_C, "%6.4lf%20.12lf\n", T, C);
        fprintf(fp_m, "%6.4lf%20.12lf\n", T, per_mag_moment);
        fprintf(fp_chi, "%6.4lf%20.12lf\n", T, chi);

        T -= 0.03;
    }

    fp_conf = fopen("conf.txt", "w");
    for(i=0; i<LENGTH; i++)
    {
        for(j=0; j<LENGTH; j++)
            fprintf(fp_conf, "%s", field[i*LENGTH + j] == 1 ? "+" : "-");
        fprintf(fp_conf, "\n");
    }
    fclose(fp_conf);

    time2 = clock();
    printf("%14.6lf s\n", (double)(time2 - time1) / CLOCKS_PER_SEC);
    fprintf(fp, "%14.6lf s\n", (double)(time2 - time1) / CLOCKS_PER_SEC);
    fclose(fp);
    fclose(fp_E);
    fclose(fp_C);
    fclose(fp_m);
    fclose(fp_chi);
    return 0;
}
void runstep(void)
{
    delta_E = 0;
    x = GetRandomInt(LENGTH);
    y = GetRandomInt(LENGTH);

    if(x != 0)
        delta_E += field[y*LENGTH + x] * field[y*LENGTH + x - 1] * 2;
    else
        delta_E += field[y*LENGTH + x] * field[y*LENGTH + LENGTH - 1] * 2;
    if(x != LENGTH)
        delta_E += field[y*LENGTH + x] * field[y*LENGTH + x + 1] * 2;
    else
        delta_E += field[y*LENGTH + x] * field[y*LENGTH] * 2;

    if(y != 0)
        delta_E += field[y*LENGTH + x] * field[y*LENGTH - LENGTH + x] * 2;
    else
        delta_E += field[y*LENGTH + x] * field[(LENGTH - 1) * LENGTH + x] * 2;
    if(y != LENGTH)
        delta_E += field[y*LENGTH + x] * field[(y+1)*LENGTH + x] * 2;
    else
        delta_E += field[y*LENGTH + x] * field[x] * 2;

    if(delta_E <= 0)
        field[y*LENGTH + x] *= -1;
    else if(GetRandomDoub() < exp(-delta_E / T))
        field[y*LENGTH + x] *= -1;
}
double GetTotalEner(void)
{
    int i, j;
    double energy = 0;
    for(i=0; i<LENGTH; i++)
    {
        for(j=0; j<LENGTH; j++)
        {
            if(i != LENGTH - 1)
                energy += - field[i*LENGTH + j] * field[(i+1) * LENGTH + j];
            else
                energy += - field[i*LENGTH + j] * field[j];

            if(j != LENGTH - 1)
                energy += - field[i*LENGTH + j] * field[i*LENGTH + j + 1];
            else
                energy += - field[i*LENGTH + j] * field[i*LENGTH];
        }
    }
    return energy;
}
double GetTotalMagMom(void)
{
    int i, j;
    double mag_mom = 0;
    for(i=0; i<LENGTH * LENGTH; i++)
        mag_mom += field[i];
    return abs(mag_mom);
}
int GetRandomInt(int range)
{
    int rndnmbr = RAND_MAX;
    while(rndnmbr > prmax)
        rndnmbr = rand();
    return rndnmbr % range;
}
double GetRandomDoub(void)
{
    return (double) rand() / RAND_MAX;
}
