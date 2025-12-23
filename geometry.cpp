#include "geometry.h"
#include <iostream>   // для cout
#include <cmath>      // для sqrt

//                 Vec

// Инициализация координат через список инициализации
Vec::Vec(double X,double Y,double Z) : x(X), y(Y), z(Z) {}

// Сложение векторов: (x,y,z) + (vx,vy,vz)
Vec Vec::add(const Vec& v) const {
    return Vec(x + v.x, y + v.y, z + v.z);
}

// Вычитание векторов: (x,y,z) - (vx,vy,vz)
Vec Vec::sub(const Vec& v) const {
    return Vec(x - v.x, y - v.y, z - v.z);
}

// Умножение на число k: (x*k, y*k, z*k)
Vec Vec::mul(double k) const {
    return Vec(x * k, y * k, z * k);
}

// Скалярное произведение: x*vx + y*vy + z*vz
double Vec::dot(const Vec& v) const {
    return x*v.x + y*v.y + z*v.z;
}

// Векторное произведение 
Vec Vec::cross(const Vec& v) const {
    return Vec(
        y*v.z - z*v.y,
        z*v.x - x*v.z,
        x*v.y - y*v.x
    );
}

// Квадрат длины |v|^2 = v•v
double Vec::norm2() const {
    return dot(*this);
}

// Длина |v| = sqrt(|v|^2)
double Vec::norm() const {
    return std::sqrt(norm2());
}

// Печать координат
void Vec::print() const {
    std::cout << "(" << x << "," << y << "," << z << ")";
}

//                  Point

// Point наследует Vec, поэтому просто вызываем конструктор Vec
Point::Point(double X,double Y,double Z) : Vec(X,Y,Z) {}

Line::Line() : ownA(0,0,0), ownB(1,0,0), a(&ownA), b(&ownB) {}

Line::Line(const Point& A, const Point& B)
    : ownA(A), ownB(B), a(&ownA), b(&ownB) {}

Line::Line(const Point* A, const Point* B)
    : ownA(), ownB(), a(A), b(B) {}

Vec Line::p() const { return *a; }
Vec Line::q() const { return *b; }
Vec Line::dir() const { return q().sub(p()); }

//                  Cramer
bool Cramer::solve(double a11,double a12,double a21,double a22,
                      double b1,double b2,double &x,double &y) {
    double D = a11*a22 - a12*a21;
    if (std::fabs(D) < 1e-12) return false;
    x = (b1*a22 - a12*b2) / D;
    y = (a11*b2 - b1*a21) / D;
    return true;
}

//                 LineChecks
double LineChecks::eps() { return 1e-9; }
bool LineChecks::isZero(double v) { return std::fabs(v) < eps(); }

bool LineChecks::parallel(const Line& l1, const Line& l2) {
    Vec u = l1.dir();
    Vec v = l2.dir();
    return u.cross(v).norm2() < eps()*eps();
}

bool LineChecks::coincident(const Line& l1, const Line& l2) {
    Vec u = l1.dir();
    Vec w = l2.p().sub(l1.p());
    return isZero(w.cross(u).norm());
}

bool LineChecks::coplanar(const Line& l1, const Line& l2) {
    Vec u = l1.dir(), v = l2.dir();
    Vec w = l2.p().sub(l1.p());
    return isZero(w.dot(u.cross(v)));
}

//                  Parallel distance 
double LineParallelOps::distance(const Line& l1, const Line& l2) {
    Vec u = l1.dir();
    Vec w = l2.p().sub(l1.p());
    return w.cross(u).norm() / u.norm();
}

static double coord(const Vec& v, int k) {
    if (k==0) return v.x; if (k==1) return v.y; return v.z;
}

static void bestPair(const Vec& n, int &i, int &j) {
    double ax=std::fabs(n.x), ay=std::fabs(n.y), az=std::fabs(n.z);
    if (ax>=ay && ax>=az) { i=1; j=2; return; }
    if (ay>=ax && ay>=az) { i=0; j=2; return; }
    i=0; j=1;
}

bool LineIntersection::compute(const Line& l1, const Line& l2, Vec& P) {
    Vec u = l1.dir(), v = l2.dir(), w = l2.p().sub(l1.p());
    Vec n = u.cross(v); int i,j; bestPair(n,i,j);

    double a11=coord(u,i), a12=-coord(v,i);
    double a21=coord(u,j), a22=-coord(v,j);
    double b1 =coord(w,i), b2 =coord(w,j);

    double t,s;
    if (!Cramer::solve(a11,a12,a21,a22,b1,b2,t,s)) return false;
    P = l1.p().add(u.mul(t));
    return true;
}

//                          Skew
bool LineSkewOps::closestPoints(const Line& l1, const Line& l2, Vec& A, Vec& B) {
    Vec u=l1.dir(), v=l2.dir(), w=l2.p().sub(l1.p());
    double a11=u.dot(u), a12=-u.dot(v), a21=u.dot(v), a22=-v.dot(v);
    double b1=u.dot(w),  b2=v.dot(w);

    double t,s;
    if (!Cramer::solve(a11,a12,a21,a22,b1,b2,t,s)) return false;
    A = l1.p().add(u.mul(t));
    B = l2.p().add(v.mul(s));
    return true;
}

void LineSkewOps::printPerpendicular(const Vec& A, const Vec& B) {
    Vec d = B.sub(A);
    std::cout << "Уравнение общего перпендикуляра:\n";
    std::cout << "x = " << A.x << " + t*" << d.x << "\n";
    std::cout << "y = " << A.y << " + t*" << d.y << "\n";
    std::cout << "z = " << A.z << " + t*" << d.z << "\n";
}

void LineTask::readPoints(Point& A, Point& B, Point& C, Point& D) {
    std::cout << "Введите A(x y z): "; std::cin >> A.x >> A.y >> A.z;
    std::cout << "Введите B(x y z): "; std::cin >> B.x >> B.y >> B.z;
    std::cout << "Введите C(x y z): "; std::cin >> C.x >> C.y >> C.z;
    std::cout << "Введите D(x y z): "; std::cin >> D.x >> D.y >> D.z;
}

void LineTask::printCoincident() { std::cout << "Совпадают (dist=0)\n"; }

void LineTask::printParallel(double dist) {
    std::cout << "Параллельны\nРасстояние = " << dist << "\n";
}

void LineTask::printIntersect(const Vec& P) {
    std::cout << "Пересекаются\nТочка пересечения = "; P.print(); std::cout << "\n";
}

void LineTask::printSkew(const Vec& A, const Vec& B) {
    std::cout << "Скрещиваются\nМинимальное расстояние = " << B.sub(A).norm() << "\n";
    std::cout << "Ближайшие точки:\nA = "; A.print(); std::cout << "\nB="; B.print(); std::cout << "\n";
    LineSkewOps::printPerpendicular(A,B);
}

int LineTask::run() {
    Point A,B,C,D; readPoints(A,B,C,D);
    Line l1(A,B); Line l2(&C,&D);

    if (LineChecks::parallel(l1,l2)) {
        if (LineChecks::coincident(l1,l2)) { printCoincident(); return 0; }
        printParallel(LineParallelOps::distance(l1,l2)); return 0;
    }
    if (LineChecks::coplanar(l1,l2)) {
        Vec P; if (LineIntersection::compute(l1,l2,P)) printIntersect(P);
        else std::cout << "Не удалось вычислить\n";
        return 0;
    }
    Vec A0,B0; if (LineSkewOps::closestPoints(l1,l2,A0,B0)) printSkew(A0,B0);
    else std::cout << "Не удалось вычислить\n";
    return 0;
}
