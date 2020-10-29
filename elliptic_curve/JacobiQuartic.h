#ifndef JacobiQuartic_h
#define JacobiQuartic_h

#include <stdio.h>
#include <openssl/bn.h>

struct MyPoint {
    BIGNUM * X;
    BIGNUM * Y;
    BIGNUM * Z;
};

struct JacobiCurve {
    // модуль, задан в parameters.h
    BIGNUM * p;
    
    // тета, найденное с помощью вольфрама
    BIGNUM * theta;
    
    // пареметры для кривой в Форме Якоби в проективной форме
    BIGNUM * eJac;
    BIGNUM * dJac;
    
    // параметры для кривой в форме Вейерштрасса в аффинных координатах, заданы в parameters.h
    BIGNUM * aCanon;
    BIGNUM * bCanon;
    BIGNUM * xJac;
    BIGNUM * yJac;
    BIGNUM * zJac;
    
    // точка на кривой в форме Вейерштрасса в проективной форме, заданы в parameters.h
    BIGNUM * xCanon;
    BIGNUM * yCanon;
    
    // порядок группы точек эллиптической кривой
    BIGNUM * q;
    
};

void initJacobiCurve(struct JacobiCurve * jacobiCurve);

// создание нейтрального элемента
void createNeutralElement( struct MyPoint * myPoint);

// создание точки
void createCustomPoint(struct MyPoint * myPoint,const char * X,const char * Y,const char * Z);

// сложение двух точек
void twoPointsAddition(struct MyPoint * results,const struct MyPoint * point1,const struct MyPoint * point2,const struct JacobiCurve * jacobiQuarticProjective);

// вывод точки в проективных координатах
void printPointProjective(const struct MyPoint * myPoint);

// вывод просто BIGNUM
void printBignum(const BIGNUM * bn);

// освобождение точки
void freeMyPoint(struct MyPoint * myPoint);

// освобождает кривую Якоби
void freeJacobiQuartic(struct JacobiCurve * jacCurve);

// проверяет, лежит ли точка на кривой; возвращает 0, если лежит и 1, -1 если не лежит
int checkPointIsOnCurve(const struct MyPoint * point,const struct JacobiCurve * jacCurve);

// сложение двух точек с использованием алгоритма Лесенки Монтгомери
void montgomeryLadder(struct MyPoint * result,const struct MyPoint * myPoint,const BIGNUM * pow,const struct JacobiCurve * jacobiCurve);

// перевод в аффинные координаты
void toAffineCoordinates (struct MyPoint *res, const struct MyPoint * myPoint, const struct JacobiCurve * jacobiCurve);

// создание базовой точки
void createBasePoint(struct MyPoint * myPoint,const struct JacobiCurve * jacobiCurve);

// вывод точки в аффинных координатах
void printPointAffine(const struct MyPoint * mypoint);

// сравнение двух точек на равенство
int twoPointsComparison(const struct MyPoint * myPoint1, const struct MyPoint * myPoint2, const struct JacobiCurve * jacobiCurve);

// to = - from
void createNegativePoint(struct MyPoint * to,const struct MyPoint * from);

#endif /* JacobiQuartic_h */
