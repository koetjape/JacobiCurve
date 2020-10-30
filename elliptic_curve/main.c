#include <stdio.h>
#include <openssl/bn.h>
#include "JacobiQuartic.h"
#include "parameters.h"
#define maxrand "10000000000000000000000000000000000000000000000000000000000000000000000"

void printResultPointOnCurve(int result){
    if (result == 0){
        printf("Точка лежит на кривой\n");
    }
    else {
        printf("Точка не лежит на кривой\n");
    }
}

void printResultPointsComparison(int result){
    if (result == 0){
        printf("Точки совпадают\n");
    }
    else {
        printf("Точки не совпадают\n");
    }
}

int main(int argc, const char * argv[]) {
    
    // инициализируем структуру jacobiCurve и наполняем ее значениями
    struct JacobiCurve jacobiCurve = {};
    initJacobiCurve(&jacobiCurve);
    
    // точка P
    struct MyPoint P = {};
    createBasePoint(&P, &jacobiCurve);
    
    printf("############################################\n\n");
    printf("e = ");
    printBignum(jacobiCurve.eJac);
    printf("d = ");
    printBignum(jacobiCurve.dJac);
    printf("p = ");
    printBignum(jacobiCurve.p);
    printf("Координаты порождающей точки в проективных координатах для кривой в форме квадрики Якоби: \n");
    printPointProjective(&P);
    printf("\n############################################\n\n");
    
    
    BIGNUM * var1 = BN_new();
    struct MyPoint var2 = {};
    struct MyPoint var3 = {};
    BIGNUM * var4 = BN_new();
    struct MyPoint var5 = {};
    // задаем верхнюю границу для генерации случайных чисел
    // генерируемые числа будут в диапазон 0 <= сгенерированное число < maxrand
    BIGNUM * high = BN_new();
    BN_dec2bn(&high, maxrand);

    // Test1
    printf("##### TEST 1 #####\n\n");
    printf("Находим Q = [k]P и проверяем, лежит точка на кривой или нет\n\n");
    // генерируем случайное число и помещаем его в var1
    BN_rand_range(var1, high);
    // выводим его:
    printf("Наша сгенерированное к:\n");
    printf("k   =   ");
    printBignum(var1);
    // подсчитываем [k]P
    montgomeryLadder(&var2, &P, var1, &jacobiCurve);
    // выводим результаты
    printPointProjective(&var2);
    toAffineCoordinates(&var3, &var2, &jacobiCurve);
    printPointAffine(&var3);
    printResultPointOnCurve(checkPointIsOnCurve(&var2, &jacobiCurve));
    
    
    
    // Test2
    printf("\n##### TEST 2 #####\n\n");
    printf("Проверяем, что [q]P = нейтральный элемент\n\n");
    printf("Точка [q]P имеет следующие координаты:\n");
    montgomeryLadder(&var2, &P, jacobiCurve.q, &jacobiCurve);
    printPointProjective(&var2);
    toAffineCoordinates(&var3, &var2, &jacobiCurve);
    printPointAffine(&var3);
    createNeutralElement(&var3);
    printResultPointsComparison(twoPointsComparison(&var3, &var2, &jacobiCurve));
    
    
    // Test3
    printf("\n##### TEST 3 #####\n\n");
    printf("Проверяем, что [q+1]P = P\n\n");
    BN_dec2bn(&var1, "1");
    BN_add(var1, var1, jacobiCurve.q);
    montgomeryLadder(&var2, &P, var1, &jacobiCurve);
    printf("Точка [q+1]P имеет следующие координаты:\n");
    printPointProjective(&var2);
    toAffineCoordinates(&var3, &var2, &jacobiCurve);
    printPointAffine(&var3);
    printf("Точка Р имеет следующие координаты:\n");
    printPointProjective(&P);
    toAffineCoordinates(&var3, &P, &jacobiCurve);
    printPointAffine(&var3);
    printResultPointsComparison(twoPointsComparison(&var2, &P, &jacobiCurve));
    
    
    // Test4
    printf("\n##### TEST 4 #####\n\n");
    printf("Проверим, что [q-1]P = -P\n\n");
    BN_dec2bn(&var1, "1");
    // var1 = q - 1
    BN_sub(var1, jacobiCurve.q, var1);
    montgomeryLadder(&var2, &P, var1, &jacobiCurve);
    printf("Точка [q-1]P имеет следующие координаты:\n");
    printPointProjective(&var2);
    toAffineCoordinates(&var3, &var2, &jacobiCurve);
    printPointAffine(&var3);
    printf("Точка -P имеет следующие координаты:\n");
    createNegativePoint(&var5, &P);
    printPointProjective(&var5);
    toAffineCoordinates(&var3, &var5, &jacobiCurve);
    printPointAffine(&var3);
    printResultPointsComparison(twoPointsComparison(&var5, &var2, &jacobiCurve));
    
    printf("\n##### TEST 5 #####\n\n");
    printf("Проверим, что [k1]P + [k2]P = [k1+k2]P\n\n");
    // var1 = k1
    BN_rand_range(var1, high);
    // var2 = [k1]P
    montgomeryLadder(&var2, &P, var1, &jacobiCurve);
    // var4 = [k2]
    BN_rand_range(var4, high);
    // var3 = [k2]P
    montgomeryLadder(&var3, &P, var4, &jacobiCurve);
    // var2 = [k1]P + [k2]P
    twoPointsAddition(&var2, &var2, &var3, &jacobiCurve);
    printf("Точка [k1]P + [k2]P имеет следующие координаты:\n");
    printPointProjective(&var2);
    toAffineCoordinates(&var3, &var2, &jacobiCurve);
    printPointAffine(&var3);
    
    // var1 = k1 + k2
    BN_add(var1,var1,var4);
    // var3 = [k1+k2]P
    montgomeryLadder(&var3, &P, var1, &jacobiCurve);
    printf("Точка [k1+k2]P имеет следующие координаты:\n");
    printPointProjective(&var3);
    toAffineCoordinates(&var5, &var3, &jacobiCurve);
    printPointAffine(&var5);
    
    // непосредственно проверяем, что [k1]P + [k2]P = [k1 + k2]P
    printResultPointsComparison(twoPointsComparison(&var2, &var3, &jacobiCurve));
    
    
    BN_free(var1);
    freeMyPoint(&var2);
    freeMyPoint(&var3);
    BN_free(var4);
    freeMyPoint(&var5);
    freeMyPoint(&P);
    freeJacobiQuartic(&jacobiCurve);
    
    
    return 0;
}
