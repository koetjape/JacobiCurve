//
//  JacobiQuartic.c
//  elliptic_curve
//
//  Created by Андрей Фомин on 23.10.2020.
//

#include <stdio.h>
#include "JacobiQuartic.h"
#include "parameters.h"

void createNeutralElement(struct MyPoint * myPoint){
    BN_dec2bn(&myPoint->X, "0");
    BN_dec2bn(&myPoint->Y, "1");
    BN_dec2bn(&myPoint->Z, "1");
}

void createCustomPoint(struct MyPoint * myPoint,const char * X,const char * Y,const char * Z){
    BN_dec2bn(&myPoint->X,X);
    BN_dec2bn(&myPoint->Y,Y);
    BN_dec2bn(&myPoint->Z,Z);
}

void printPointProjective(const struct MyPoint * myPoint){
    printf("Проективные координаты:\nX   = %s\nY   = %s\nZ   = %s\n", BN_bn2dec(myPoint->X), BN_bn2dec(myPoint->Y), BN_bn2dec(myPoint->Z));
    
}

// печатает число, хранящееся в BIGNUM в дестичной системе
void printBignum(const BIGNUM * bn){
    printf("%s\n", BN_bn2dec(bn));
}

// освобождает память точки
void freeMyPoint(struct MyPoint * myPoint){
    BN_free(myPoint->X);
    BN_free(myPoint->Y);
    BN_free(myPoint->Z);
}

// производит сложение двух точек poin1 и point2 на эллиптической кривой jacobiQuarticProjective, результат помещается в result
void twoPointsAddition(struct MyPoint * result,const struct MyPoint * point1,const struct MyPoint * point2,const struct JacobiCurve * jacobiQuartic){
    // T1 = X1; T2 = Y1; T3 = Z1; T4 = X2; T5 = Y2; T6 = Z2
    BN_CTX * bnCtx = BN_CTX_new();
    BIGNUM * T1 = BN_dup(point1->X);
    BIGNUM * T2 = BN_dup(point1->Y);
    BIGNUM * T3 = BN_dup(point1->Z);
    BIGNUM * T4 = BN_dup(point2->X);
    BIGNUM * T5 = BN_dup(point2->Y);
    BIGNUM * T6 = BN_dup(point2->Z);
    BIGNUM * T7 = BN_new();
    BIGNUM * T8 = BN_new();
    
//    printBignum(T1);
//    printBignum(T2);
//    printBignum(T3);
//    printBignum(T4);
//    printBignum(T5);
//    printBignum(T6);

    // T7 = T1 * T3
    BN_mod_mul(T7, T1, T3, jacobiQuartic->p, bnCtx);
    // T7 = T2 + T7
    BN_mod_add(T7,T2,T7,jacobiQuartic->p,bnCtx);
    // T8 = T4 * T6
    BN_mod_mul(T8,T4,T6,jacobiQuartic->p,bnCtx);
    // T8 = T5 + T8
    BN_mod_add(T8,T5,T8, jacobiQuartic->p,bnCtx);
    // T2 = T2 * T5
    BN_mod_mul(T2, T2,T5,jacobiQuartic->p,bnCtx);
    // T7 = T7 * T8
    BN_mod_mul(T7,T7,T8,jacobiQuartic->p,bnCtx);
    // T7 = T7 - T2
    BN_mod_sub(T7,T7,T2,jacobiQuartic->p, bnCtx);
    // T5 = T1 * T4
    BN_mod_mul(T5,T1,T4,jacobiQuartic->p, bnCtx);
    // T1 = T1 + T3
    BN_mod_add(T1,T1,T3,jacobiQuartic->p, bnCtx);
    // T8 = T3 * T6
    BN_mod_mul(T8,T3,T6,jacobiQuartic->p,bnCtx);
    // T4 = T4 + T6
    BN_mod_add(T4,T4,T6, jacobiQuartic->p,bnCtx);
    // T6 = T5 * T8
    BN_mod_mul(T6, T5, T8, jacobiQuartic->p, bnCtx);
    // T7 = T7 - T6
    BN_mod_sub(T7,T7,T6,jacobiQuartic->p,bnCtx);
    // T1 = T1 * T4
    BN_mod_mul(T1,T1,T4,jacobiQuartic->p,bnCtx);
    // T1 = T1 - T5
    BN_mod_sub(T1,T1,T5,jacobiQuartic->p, bnCtx);
    // T1 = T1 - T8
    BN_mod_sub(T1,T1,T8, jacobiQuartic->p, bnCtx);
    // T3 = T1 * T1
    BN_mod_mul(T3,T1,T1,jacobiQuartic->p,bnCtx);
    // T6 = T6 + T6
    BN_mod_add(T6,T6,T6,jacobiQuartic->p,bnCtx);
    // T3 = T3 - T6
    BN_mod_sub(T3,T3,T6,jacobiQuartic->p,bnCtx);
    // T4 = e * T6
    BN_mod_mul(T4, jacobiQuartic->eJac,T6,jacobiQuartic->p, bnCtx);
    // T3 = T3 * T4
    BN_mod_mul(T3,T3,T4,jacobiQuartic->p,bnCtx);
    // T4 = d * T6
    BN_mod_mul(T4,jacobiQuartic->dJac,T6,jacobiQuartic->p,bnCtx);
    // T2 = T2 - T4
    BN_mod_sub(T2,T2,T4, jacobiQuartic->p, bnCtx);
    // T4 = T8 * T8
    BN_mod_mul(T4,T8,T8,jacobiQuartic->p,bnCtx);
    // T8 = T5 * T5
    BN_mod_mul(T8,T5,T5,jacobiQuartic->p,bnCtx);
    // T8 = e * T8
    BN_mod_mul(T8,jacobiQuartic->eJac,T8,jacobiQuartic->p,bnCtx);
    // T5 = T4 + T8
    BN_mod_add(T5,T4,T8,jacobiQuartic->p,bnCtx);
    // T2 = T2 * T5
    BN_mod_mul(T2,T2,T5,jacobiQuartic->p,bnCtx);
    // T2 = T2 + T3
    BN_mod_add(T2,T2,T3,jacobiQuartic->p,bnCtx);
    // T5 = T4 - T8
    BN_mod_sub(T5,T4,T8,jacobiQuartic->p,bnCtx);
    result->X = BN_dup(T7);
    result->Y = BN_dup(T2);
    result->Z = BN_dup(T5);
    BN_free(T1);
    BN_free(T2);
    BN_free(T3);
    BN_free(T4);
    BN_free(T5);
    BN_free(T6);
    BN_free(T7);
    BN_free(T8);
    BN_CTX_free(bnCtx);
    
}

void freeJacobiQuartic(struct JacobiCurve * jacobiCurve){
    BN_free(jacobiCurve->dJac);
    BN_free(jacobiCurve->eJac);
    BN_free(jacobiCurve->p);
    BN_free(jacobiCurve->xJac);
    BN_free(jacobiCurve->yJac);
    BN_free(jacobiCurve->zJac);
    BN_free(jacobiCurve->aCanon);
    BN_free(jacobiCurve->bCanon);
    BN_free(jacobiCurve->xCanon);
    BN_free(jacobiCurve->yCanon);
    BN_free(jacobiCurve->q);
}

int checkPointIsOnCurve(const struct MyPoint * point,const  struct JacobiCurve * jacobiCurve) {
    BIGNUM * temp1 = BN_new();
    BIGNUM * temp2 = BN_new();
    BN_CTX * bnCtx = BN_CTX_new();
    BIGNUM * four = BN_new();
    BN_dec2bn(&four, "4");
    
    // temp1 = X * Z
    BN_mod_mul(temp1, point->X, point->Z, jacobiCurve->p, bnCtx);
    // temp1 = X^2 * Z^2
    BN_mod_sqr(temp1, temp1, jacobiCurve->p, bnCtx);
    // temp1 = d * X^2 * Z^2
    BN_mod_mul(temp1, temp1, jacobiCurve->dJac, jacobiCurve->p, bnCtx);
    // temp1 = 2 * d * X^2 * Z^2 = d * X^2 * Z^2 + d * X^2 * Z^2
    BN_mod_add(temp1, temp1, temp1, jacobiCurve->p, bnCtx);
    // temp2 = X^4
    BN_mod_exp(temp2, point->X,four, jacobiCurve->p, bnCtx);
    // temp2 = e * X^4
    BN_mod_mul(temp2, temp2, jacobiCurve->eJac, jacobiCurve->p, bnCtx);
    // temp2 = e * X^4 - 2 * d * X^2 * Z^2
    BN_mod_sub(temp2, temp2, temp1, jacobiCurve->p, bnCtx);
    // temp1 = Z^4
    BN_mod_exp(temp1, point->Z,four, jacobiCurve->p, bnCtx);
    // temp2 = e * X^4 - 2 * d * X^2 * Z^2 + Z^4
    BN_mod_add(temp2, temp2, temp1, jacobiCurve->p, bnCtx);
    // temp1 = Y^2
    BN_mod_sqr(temp1, point->Y, jacobiCurve->p, bnCtx);
    
    int result = BN_cmp(temp1, temp2);
    
    BN_free(temp1);
    BN_free(temp2);
    BN_free(four);
    BN_CTX_free(bnCtx);
    
    return result;
}

void montgomeryLadder(struct MyPoint * Q,const struct MyPoint * P,const BIGNUM * k,const struct JacobiCurve * jacobiCurve){
    int numberOfBits = BN_num_bits(k);
    struct MyPoint R;
    // R = P
    R.X = BN_dup(P->X);
    R.Y = BN_dup(P->Y);
    R.Z = BN_dup(P->Z);
    // Q = нейтральный элемент
    createNeutralElement(Q);
    for (int i = numberOfBits; i > -1 ;i--){
        if (BN_is_bit_set(k, i)){
            twoPointsAddition(Q, Q, &R, jacobiCurve);
            twoPointsAddition(&R, &R, &R, jacobiCurve);

        }
        else {
            twoPointsAddition(&R, &R, Q, jacobiCurve);
            twoPointsAddition(Q, Q, Q, jacobiCurve);
        }
    }
    
    freeMyPoint(&R);
    
}

void initJacobiCurve(struct JacobiCurve * jacobiCurve){
    BN_dec2bn(&jacobiCurve->p, pPar);
    BN_dec2bn(&jacobiCurve->aCanon, aPar);
    BN_dec2bn(&jacobiCurve->bCanon, bPar);
    BN_dec2bn(&jacobiCurve->xCanon, xPar);
    BN_dec2bn(&jacobiCurve->yCanon, yPar);
    BN_dec2bn(&jacobiCurve->theta, thetaPar);
    BN_dec2bn(&jacobiCurve->q, qPar);
    
    BIGNUM * temp1 = BN_new();
    BIGNUM * temp2 = BN_new();
    BN_CTX * bnCtx = BN_CTX_new();
    
    // находим d (параметр кривой в форме квартики Якоби в проективных координатах)
    // temp1 = 3
    BN_dec2bn(&temp1, "3");
    // temp1 = 3 * theta
    BN_mod_mul(temp1, temp1, jacobiCurve->theta, jacobiCurve->p, bnCtx);
    // temp2 = 4
    BN_dec2bn(&temp2, "4");
    // temp2 = 4^(-1)
    BN_mod_inverse(temp2, temp2, jacobiCurve->p, bnCtx);
    // temp1 = 3 * theta * 4^(-1) = d
    BN_mod_mul(temp1, temp1, temp2, jacobiCurve->p, bnCtx);
    
    jacobiCurve->dJac = BN_dup(temp1);
    // теперь превратим параметр d  в параметр e с помощью следующих преобразований
    // temp1 = 3 * theta * 4^(-1) * theta = 3 * theta^2 * 4^(-1)
    BN_mod_mul(temp1, temp1, jacobiCurve->theta, jacobiCurve->p, bnCtx);
    // temp1 = 3 * theta^2 * 4^(-1) + a
    BN_mod_add(temp1, temp1, jacobiCurve->aCanon, jacobiCurve->p, bnCtx);
    // temp1 =( 3 * theta^2 * 4^(-1) + a ) * 4^(-1) = e
    BN_mod_mul(temp1, temp1, temp2, jacobiCurve->p, bnCtx);
    // temp2 = -1
    BN_dec2bn(&temp2, "-1");
    // temp1 = -1 * (3 * theta^2 * 4^(-1) + a) * 4^(-1)
    BN_mod_mul(temp1, temp1, temp2, jacobiCurve->p, bnCtx);
    
    jacobiCurve->eJac = BN_dup(temp1);
    
    
    // переведем (x,y) в (X, Y, Z)
    // temp1 = (x - theta)
    BN_mod_sub(temp1, jacobiCurve->xCanon, jacobiCurve->theta, jacobiCurve->p, bnCtx);
    // temp2 = 2
    BN_dec2bn(&temp2, "2");
    // temp1 = 2 * (x - theta) = X
    BN_mod_mul(temp1, temp1, temp2, jacobiCurve->p, bnCtx);
    jacobiCurve->xJac = BN_dup(temp1);
    
    // temp1 = x - theta
    BN_mod_sub(temp1, jacobiCurve->xCanon, jacobiCurve->theta, jacobiCurve->p, bnCtx);
    // temp1 = (x - theta)^2
    BN_mod_sqr(temp1, temp1, jacobiCurve->p, bnCtx);
    // temp2 = x + x = 2x
    BN_mod_add(temp2, jacobiCurve->xCanon, jacobiCurve->xCanon, jacobiCurve->p, bnCtx);
    // temp2 = 2x + theta
    BN_mod_add(temp2, temp2, jacobiCurve->theta, jacobiCurve->p, bnCtx);
    // temp1 = (2x + theta) (x - theta)^2
    BN_mod_mul(temp1, temp2, temp1, jacobiCurve->p, bnCtx);
    // temp2 = y^2
    BN_mod_sqr(temp2, jacobiCurve->yCanon, jacobiCurve->p, bnCtx);
    // temp1 = (2x + theta) (x - theta)^2 - y^2 = Y
    BN_mod_sub(temp1, temp1, temp2, jacobiCurve->p, bnCtx);
    jacobiCurve->yJac = BN_dup(temp1);
    // Z = y
    jacobiCurve->zJac = BN_dup(jacobiCurve->yCanon);
    
    // порядок группы m
    BN_dec2bn(&jacobiCurve->q, qPar);
    
    BN_CTX_free(bnCtx);
    BN_free(temp1);
    BN_free(temp2);
}

// перевод точки в аффинные координаты
void toAffineCoordinates (struct MyPoint *res, const struct MyPoint *myPoint, const struct JacobiCurve *jacobiCurve)
{
    BIGNUM *x = BN_new(),
           *y = BN_new(),
           *z = BN_new();

    BN_CTX *tmp = BN_CTX_new ();

    // x = Z^(-1)
    BN_mod_inverse (x, myPoint->Z, jacobiCurve->p, tmp);
    // x = X / Z
    BN_mod_mul (x, x, myPoint->X, jacobiCurve->p, tmp);

    // y = 1 / Z
    BN_mod_inverse (y, myPoint->Z, jacobiCurve->p, tmp);
    // y = 1 / Z^2
    BN_mod_mul (y, y, y, jacobiCurve->p, tmp);
    // y = Y / Z^2
    BN_mod_mul (y, y, myPoint->Y, jacobiCurve->p, tmp);

    // z = 0
    BN_dec2bn (&z, "0");

    res->X = BN_dup(x);
    res->Y = BN_dup(y);
    res->Z = BN_dup(z);

    BN_free (x);
    BN_free (y);
    BN_free (z);
    BN_CTX_free (tmp);
}
// делаем myPoint "базовой" точкой в проективных координатах
void createBasePoint(struct MyPoint * myPoint,const struct JacobiCurve * jacobiCurve){
    myPoint->X = BN_dup(jacobiCurve->xJac);
    myPoint->Y = BN_dup(jacobiCurve->yJac);
    myPoint->Z = BN_dup(jacobiCurve->zJac);
}

void printPointAffine(const struct MyPoint * myPoint){
    printf("Аффинные координаты:\nx   = %s\ny   = %s\n", BN_bn2dec(myPoint->X), BN_bn2dec(myPoint->Y));
}

int twoPointsComparison(const struct MyPoint * myPoint1, const struct MyPoint * myPoint2,const struct JacobiCurve * jacobiCurve){
    struct MyPoint temp1 = {};
    struct MyPoint temp2 = {};
    int result;
    
    toAffineCoordinates(&temp1, myPoint1, jacobiCurve);
    toAffineCoordinates(&temp2, myPoint2, jacobiCurve);
    if (BN_cmp(temp1.X, temp2.X) == 0 && BN_cmp(temp1.Y, temp1.Y) == 0){
        result = 0;
    }
    else {
        result = -1;
    }
    return result;
}
// создаем противоположную точку
// from : (X, Y, Z), to : (-X, Y, Z)
void createNegativePoint(struct MyPoint * to, const struct MyPoint * from){
    to->X = BN_dup(from->X);
    to->Y = BN_dup(from->Y);
    to->Z = BN_dup(from->Z);
    BN_set_negative(to->X, 1);
}


