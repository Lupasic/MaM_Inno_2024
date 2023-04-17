/*
С помощью  неявной разностной схемы решить нестационарное уравнение теплопроводности для прямоугольной пластины  размером  8*3см. Начальное значение температуры пластины - 10 градусов.
Граничные условия следующие:  верхняя и верхняя половина правой границы теплоизолированы, на остальной части границы температура 400 градусов.
При выводе результатов показать динамику изменения температуры (например с помощью цветовой гаммы).
     Отчет должен содержать : текст программы, рисунок объекта с распределением температуры в момент времени 20 сек сравнение результатов расчета с результатами, полученными с помощью пакета ANSYS .
*/
#define nodePerL 1  // число узлов на ед длины
#define LENX 9    // число узлов по x
#define LENY 5     // число узлов по y
#define T0 10.0       // начальная температура
#define T1 400.0      // граничное условие право-низ
#define T2 0.0 // термоизоляция
#define Leng LENX*LENY //Длина


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <fcntl.h>

int Nm, modTime;  // nm - размер матрицы Modtime - время решения (то которое задано по заданию)
double mA[Leng*Leng], vX[Leng], vB[Leng]; // вектора матрицы
double L, H, deltaX, deltaT, startT, Trl; // длина, высота, deltax = deltay
double startT, Trl; // startT - начальная температура, остальное - граничные условия
double a = 1.0; // a  из формул
 int NoX, NoY, edge; //число узлов по x и по y из дефайна

// Гаусс [5]
 void frw_one_th () //Прямой ход
{
     int i, j, k;
long double dgE;

    for (k = 0; k < Nm; k++)
    {
        dgE = mA[Nm*k + k];
        for (j = k; j < Nm; j++)
            mA[Nm*k + j]/= dgE;
        vB[k]/= dgE;

        for (i = k + 1; i < Nm; i++)
        {
            dgE = mA[Nm*i + k];
            for (j = k; j < Nm; j++)
                mA[Nm*i + j]-= mA[Nm*k + j]*dgE;
            vB[i]-= vB[k]*dgE;
        }
    }
}

void bck_one_th () //Обратный ход
{
     int i, j;

    vX[Nm - 1] = vB[Nm - 1];

    for (i = Nm - 2; i >= 0; i--)
    {
        vX[i] = vB[i];
        for (j = i + 1; j < Nm; j++)
            vX[i]-= mA[Nm*i + j]*vX[j];
    }

}
//[5]


void all_gen() //генерация начальной температуры
{
    int i;
    for (i = 0; i < Nm; i++)
    {
        vX[i] = startT;
    }
}


void  generan_matrx()
{

/*почему такое обращение к координатам :
    Nm*(sv) - перемещаемся по уравнениям(каждый раз перемещаемся на одно уравнение ниже)
    когда i +/- 1 то просто прибавляем их
    когда j+/-1 то +/- NoX так как у нас уравнение записывется в одну линию и чтобы переместиться на слой выше/ниже мы проходим как раз 1 раз по х
    NoX*j - чтобы записывать уравнения в зависимости от уровня(смещение получается как бы)
*/
int i=0, j=0, sv=0;
    //заполнение 0ями массивы
    for (i = 0; i < Nm*Nm; i++)
        mA[i] = 0.0;
    for (i = 0; i < Nm; i++)
        vB[i] = 0.0;


    for (j = 0; j < NoY; j++)
    {
        for (i = 0; i < NoX; i++)
        {

            if((i==NoX-1 && j==NoY-1)) //Углы, не понятно как их правильно (i==0 && j==NoY-1) ||
            {
                mA[Nm*(sv) + i + (NoX*j)] = 1.0/deltaX;
                   mA[Nm*(sv) + i - 1 + (NoX*j)] = -1.0/deltaX;
               vB[sv] = 0.0;
               sv++;
               continue;
            }

            if(j==NoY-1)//условие второго рода сверху
            {
                 mA[Nm*(sv) + i + (NoX*j)] = 1.0/deltaX;
                mA[Nm*(sv) + (i-NoX) + (NoX*j)] = -1.0/deltaX;
                vB[sv] = T2;
                sv++;
                continue;
            }
            if(i==NoX-1 && j<NoY/2) //условие первого рода правая половина низ
                {
                    vB[sv] = Trl;
                    mA[Nm*(sv) + i + (NoX*j)] = 1.0/deltaX;
                    sv++;
                    continue;
                }
                if(i==NoX-1 && j>=NoY/2) //условие второго рода правая половина вверх
                    {
                         mA[Nm*(sv) + i + (NoX*j)] = 1.0/deltaX;
                        mA[Nm*(sv) + i - 1 + (NoX*j)] = -1.0/deltaX;
                        vB[sv] = T2;
                        sv++;
                        continue;
                    }
                           // Стандартный случай заполнения
                        if((i-1)>=0)
                             mA[Nm*(sv) + i - 1 + (NoX*j)] = -a/deltaX*deltaX;
                        mA[Nm*(sv) + i+(NoX*j)] = 1.0/deltaT + 4.0*a/deltaX*deltaX;
                        if((i+1)<=Nm)
                             mA[Nm*(sv) + i + 1 + (NoX*j)] = -a/deltaX*deltaX;
                        if((j*NoX+i+NoX)<=Nm)
                             mA[Nm*(sv) + (i + NoX) + (NoX*j)] = -a/deltaX*deltaX;
                         if((j*NoX+i-NoX)>=0)
                             mA[Nm*(sv) + (i - NoX) + (NoX*j)] = -a/deltaX*deltaX;
                             vB[sv] = vX[sv]/deltaT;
                             sv++;                 
      }
    }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
    int i=0, j, k, ret;
    char string[80];
    FILE *fp;
   fp = fopen( "result", "w" ); //открытие файла на запись данных
    // проверка на наличие аргумента [1]
    if (argc < 2)
    {
            printf("Нету аргумента - времени\n");
            return(1);
    }
    //[1]
    modTime = atoi(argv[1]); //задача времени
    //задаем данные из дефайна и глобальных переменных [2]
    L = LENX;
    H = LENY;
    Nm = Leng; // количество узлов всего
    NoX = LENX*nodePerL;
    NoY = LENY*nodePerL;
    deltaX = H*L/Nm; // дельта X
    deltaT = 1; // дельта t
    startT = T0;
    Trl = T1;
    //[2]

    printf("Количество узлов: %d\n", Nm);
    all_gen();  // генерация начальных данных
    for (i = 0; i < modTime/deltaT; i++)	//решение уравнений
    {
        generan_matrx();
        frw_one_th();
        bck_one_th ();
        // Печать
       //ret = Nm+NoX;
        ret=0;

            for (j = 0; j < NoY; j++)
            {
              //  ret=ret-2*NoX;
                for (k = 0; k < NoX; k++)
                {
                    sprintf(string, "%.3f ",vX[ret++]);
                    fputs(string,fp);
                }
                sprintf(string, "\n");
                fputs(string,fp);
            }
            sprintf(string, "\n\n");
            fputs(string,fp);
    }
    return(0);
}
