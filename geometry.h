#pragma once

/// @brief Vec - вектор (x,y,z) ("стрелка")
/// описывает направление
class Vec{
public:
    double x,y,z;

    Vec(double X=0, double Y=0, double Z=0);

    Vec add(const Vec& p) const; //this+p
    Vec sub(const Vec& p) const; //this-p
    Vec mul(double k) const; //this*k

    double dot(const Vec& p) const; //this*p
    Vec cross(const Vec& p) const; // perpend

    double norm() const; // |p|
    double norm2() const; // |p|^2

    void print() const; // печать координат
};

/// @brief Point - точка в пространстве
class Point:public Vec{
public:
    Point(double X=0, double Y=0, double Z=0);
};


class Line{
    Point ownA, ownB; // композиция: хранение копии точек
    const Point* a; // агрегация: ссылка на внешние точки
    const Point* b;

public:
    Line();

    Line(const Point& A, const Point& B);// композиция: копируем A,B внутрь
    Line(const Point* A, const Point* B);// агрегация: держим указатели

    Vec p() const; // первая точка
    Vec q() const; // вторая
    Vec dir() const; // направление (q - p)
};

class Cramer{
public:
    static bool solve(double a11,double a12,double a21,double a22,
                      double b1,double b2,double &x,double &y);
};

class LineChecks {
public:
    static double eps();
    static bool isZero(double v);

    static bool parallel(const Line& l1, const Line& l2);
    static bool coincident(const Line& l1, const Line& l2);
    static bool coplanar(const Line& l1, const Line& l2);
};

class LineParallelOps {
public:
    static double distance(const Line& l1, const Line& l2);
};

class LineIntersection {
public:
    static bool compute(const Line& l1, const Line& l2, Vec& P);
};

class LineSkewOps {
public:
    static bool closestPoints(const Line& l1, const Line& l2, Vec& A, Vec& B);
    static void printPerpendicular(const Vec& A, const Vec& B);
};

class LineTask {
    static void readPoints(Point& A, Point& B, Point& C, Point& D);
    static void printIntersect(const Vec& P);
    static void printParallel(double dist);
    static void printCoincident();
    static void printSkew(const Vec& A, const Vec& B);
public:
    static int run();
};
