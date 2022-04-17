#define _USE_MATH_DEFINES 
#include <cmath>
#include<iostream>
#include<vector>
#include<string>
#include<algorithm>

bool comp(double a, double b) {
	return abs(a - b) < 1e-7;
}

class Line;
class Vector;

class Point {
private:
	friend Point operator+(Point start, Vector transfer);
	friend Point Middle(Point first, Point second);
	friend Point Divide(Point first, Point second, double k);
	friend double Distance(Point first, Point second);
	friend Point Projection(Point point, Line line);
	friend double Distance(Point point, Line line);

public:
	double x;
	double y;

	Point(void);
	Point(double x_, double y_);
	Point(Vector vector);
	//параллельный перенос
	Point& operator+=(Vector transfer);
	//оператор присваивания
	Point& operator=(Point another);
	//меньше та точка, которая левее и ниже
	bool operator<(Point another) const;
	//проверка совпадения
	bool operator==(Point another) const;
	//проверка несовпадения
	bool operator!=(Point another) const;
	//поворот
	void rotate(Point center, double angle);
	//центральная симметрия
	void reflex(Point center);
	//осевая симметрия
	void reflex(Line axis);
	//гомотетия
	void scale(Point center, double coefficient);
};




class Vector {
private:
	friend double operator*(Vector first, Vector second);
	friend double operator^(Vector first, Vector second);

public:
	double x;
	double y;

	Vector(void);
	//радиус-вектор точки
	Vector(Point point);
	//вектор по координатам
	Vector(double x_, double y_);
	//вектор по двум точкам
	Vector(Point first, Point second);
	//прибавление вектора
	Vector& operator+=(Vector another);
	//вычитание вектора
	Vector& operator-=(Vector another);
	//умножение на число
	Vector& operator*=(double k);
	//деление на число
	Vector& operator/=(double k);
	//проверка равенства
	bool operator==(Vector another) const;
	//проверка неравенства
	bool operator!=(Vector another) const;
	//проверка коллинеарности
	bool operator||(Vector another) const;
	//оператор присваивания
	Vector operator=(Point point);
	//длина вектора
	double Len() const;
	//поворот
	void rotate(Point center, double angle);
	//осевая симметрия
	void reflex(Line axis);

};




class Line {
protected:
	friend Point Intersect(Line first, Line second);
	friend Line Per(Point point, Line line);
	friend Line MidPer(Point first, Point second);
	friend Line Bisector(Point A, Point O, Point B);

public:
	Vector rds; //radius-vector - радиус-вектор
	Vector drct; //direction - направляющий вектор

	Line(void);
	//прямая по точке и направляющему вектору
	Line(Point point, Vector direction);
	//прямая по радиус-вектору и направляющему вектору
	Line(Vector radius, Vector direction);
	//прямая по точке и параллельной прямой
	Line(Point point, Line parallel);
	//прямая по двум точкам
	Line(Point first, Point second);
	//прямая по угловому коэффициенту и сдвигу
	Line(double k, double b);
	//прямая по точке и угловому коэффициенту
	Line(Point point, double k);
	//проверка совпадения
	bool operator==(Line another) const;
	//проверка несовпадения
	bool operator!=(Line another) const;
	//проверка параллельности
	bool operator||(Line another) const;
	//поворот
	void rotate(Point center, double angle);
	//центральная симметрия
	void reflex(Point center);
	//осевая симметрия
	void reflex(Line axis);
	//гомотетия
	void scale(Point center, double coefficient);
};




class Shape {
protected:
	Shape() = default;

public:
	//периметр
	virtual double perimeter() const = 0;
	//площадь
	virtual double area() const = 0;
	//равенство (по множествам точек)
	virtual bool operator==(const Shape& another) const = 0;
	//неравенство (по множествам точек)
	virtual bool operator!=(const Shape& another) const = 0;
	//равенство
	virtual bool isCongruentTo(const Shape& another) const = 0;
	//подобие
	virtual bool isSimilarTo(const Shape& another) const = 0;
	//содержит ли точку
	virtual bool containsPoint(Point point) const = 0;
	//поворот
	virtual void rotate(Point center, double angle) = 0;
	//центральная симметрия
	virtual void reflex(Point center) = 0;
	//осевая симметрия
	virtual void reflex(Line axis) = 0;
	//гомотетия
	virtual void scale(Point center, double coefficient) = 0;

	virtual ~Shape() = default;
};




class Ellipse : public Shape {
protected:
	Point F1;
	Point F2;
	double focal_distance;

public:
	Ellipse(void);
	Ellipse(Point F1_, Point F2_, double focal_distance_);
	//фокусы
	std::pair<Point, Point> focuses() const;
	//директрисы
	std::pair<Line, Line> directrices() const;
	//эксцентриситет
	double eccentricity() const;
	//центр 
	Point center() const;
	//F1 - левый, F2 - правый
	void Standartization();
	//периметр
	double perimeter() const override;
	//площадь
	double area() const override;
	//равенство (по множествам точек)
	bool operator==(const Shape& another) const override;
	//неравенство (по множествам точек)
	bool operator!=(const Shape& another) const override;
	//равенство
	bool isCongruentTo(const Shape& another) const override;
	//подобие
	bool isSimilarTo(const Shape& another) const override;
	//содержит ли точку
	bool containsPoint(Point point) const override;
	//поворот
	void rotate(Point center, double angle) override;
	//центральная симметрия
	void reflex(Point center) override;
	//осевая симметрия
	void reflex(Line axis) override;
	//гомотетия
	void scale(Point center, double coefficient) override;
};




class Circle : public Ellipse {
protected:

public:
	Point O;
	double R;

	Circle(void);
	//окружность по центру и радиусу
	Circle(const Point center_, double radius_);
	//радиус
	double radius() const;
	//периметр
	double perimeter() const override;
	//площадь
	double area() const override;
	//поворот
	void rotate(Point center, double angle) override;
	//центральная симметрия
	void reflex(Point center) override;
	//осевая симметрия
	void reflex(Line axis) override;
	//гомотетия
	void scale(Point center, double coefficient) override;
};




class Polygon : public Shape {
protected:
	std::vector<Point> points;

public:
	Polygon(void);
	Polygon(const std::vector<Point> points_);
	explicit Polygon(std::initializer_list<Point> list);
	//отношение соседних сторон
	double ratio_of_side(size_t i) const;
	//синус угла
	double angle(size_t i) const;
	//количество вершин
	int verticesCount() const;
	//вектор вершин
	const std::vector<Point>& getVertices() const;
	//выпуклость
	bool isConvex() const;
	//периметр
	double perimeter() const override;
	//площадь
	double area() const override;
	//равенство (по множествам точек)
	bool operator==(const Shape& another) const override;
	//неравенство (по множествам точек)
	bool operator!=(const Shape& another) const override;
	//равенство
	bool isCongruentTo(const Shape& another) const override;
	//подобие
	bool isSimilarTo(const Shape& another) const override;
	//содержит ли точку
	bool containsPoint(Point point) const override;
	//поворот
	void rotate(Point center, double angle) override;
	//центральная симметрия
	void reflex(Point center) override;
	//осевая симметрия
	void reflex(Line axis) override;
	//гомотетия
	void scale(Point center, double coefficient) override;
};




class Rectangle : public Polygon {
protected:
	Point first;
	Point second;
	double side_ratio;

public:
	Rectangle(void);
	Rectangle(Point first_, Point second_, double side_ratio_);
	//центр
	Point center() const;
	//диагонали
	std::pair<Line, Line> diagonals() const;
	//первая вершина - самая маленькая, далее в порядке обхода по часовой стрелке
	void Standartization();
	//поворот
	void rotate(Point center, double angle) override;
	//центральная симметрия
	void reflex(Point center) override;
	//осевая симметрия
	void reflex(Line axis) override;
	//гомотетия
	void scale(Point center, double coefficient) override;
};




class Square : public Rectangle {
protected:

public:
	Square(void);
	Square(Point first_, Point second_);
	//описанная окружность
	Circle circumscribedCircle() const;
	//вписанная окружность
	Circle inscribedCircle() const;
};




class Triangle : public Polygon {
protected:
	Point A;
	Point B;
	Point C;

public:
	Triangle(void);
	Triangle(Point A_, Point B_, Point C_);
	//описанная окружность
	Circle circumscribedCircle() const;
	//вписанная окружность
	Circle inscribedCircle() const;
	//центр описанной окружности
	Point circumcenter() const;
	//инцентр
	Point incenter() const;
	//центроид
	Point centroid() const;
	//ортоцентр
	Point orthocenter() const;
	//прямая Эйлера
	Line EulerLine() const;
	//окружность Эйлера/9-ти точек
	Circle ninePointsCircle() const;
	//первая вершина самая маленькая, далее в порядке обхода по часовой стрелке
	void Standartization();
	//периметр
	void rotate(Point center, double angle) override;
	//центральная симметрия
	void reflex(Point center) override;
	//осевая симметрия
	void reflex(Line axis) override;
	//гомотетия
	void scale(Point center, double coefficient) override;
};























Vector::Vector(void) = default;

//радиус-вектор точки
Vector::Vector(Point point) : x(point.x), y(point.y) {}

//вектор по координатам
Vector::Vector(double x_, double y_) : x(x_), y(y_) {}

//вектор по двум точкам
Vector::Vector(Point first, Point second) : x(second.x - first.x), y(second.y - first.y) {}

//прибавление вектора
Vector& Vector::operator+=(Vector another) {
	x += another.x;
	y += another.y;
	return *this;
}
//вычитание вектора
Vector& Vector::operator-=(Vector another) {
	x -= another.x;
	y -= another.y;
	return *this;
}
//умножение на число
Vector& Vector::operator*=(double k) {
	x *= k;
	y *= k;
	return *this;
}
//деление на число
Vector& Vector::operator/=(double k) {
	x /= k;
	y /= k;
	return *this;
}
//проверка равенства
bool Vector::operator==(Vector another) const {
	return comp(x, another.x) && comp(y, another.y);
}
//проверка неравенства
bool Vector::operator!=(Vector another) const {
	return !(*this == another);
}
//проверка коллинеарности
bool Vector::operator||(Vector another) const {
	return comp((*this ^ another), 0);
}
//оператор присваивания
Vector Vector::operator=(Point point) {
	x = point.x;
	y = point.y;
	return *this;
}
//длина вектора
double Vector::Len() const {
	double ans = std::sqrt(x * x + y * y);
	return ans;
}
//поворот
void Vector::rotate(Point center, double angle) {
	double s = std::sin(M_PI * angle / 180);
	double c = std::cos(M_PI * angle / 180);
	double x_ = x;
	double y_ = y;
	x = center.x + c * (x_ - center.x) - s * (y_ - center.y);
	y = center.y + s * (x_ - center.x) + c * (y_ - center.y);
}
//осевая симметрия
void Vector::reflex(Line axis) {
	Point begin = Point(0, 0);
	Point end = Point(x, y);
	begin.reflex(axis);
	end.reflex(axis);
	*this = Vector(begin, end);
}
//сумма векторов
Vector operator+(Vector first, Vector second) {
	Vector copy = first;
	copy += second;
	return copy;
}
//разность векторов
Vector operator-(Vector first, Vector second) {
	Vector copy = first;
	copy -= second;
	return copy;
}
//умножение вектора на число
Vector operator*(Vector vector, double k) {
	Vector copy = vector;
	copy *= k;
	return copy;
}
//скалярное произведение
double operator*(Vector first, Vector second) {
	return first.x * second.x + first.y * second.y;
}
//векторное произведение
double operator^(Vector first, Vector second) {
	return first.x * second.y - second.x * first.y;
}
//деление вектора на число
Vector operator/(Vector vector, double k) {
	Vector copy = vector;
	copy /= k;
	return copy;
}








Point::Point(void) = default;

Point::Point(double x_, double y_) : x(x_), y(y_) {}

Point::Point(Vector vector) : x(vector.x), y(vector.y) {}

//параллельный перенос
Point& Point::operator+=(Vector transfer) {
	x += transfer.x;
	y += transfer.y;
	return *this;
}
//оператор присваивания
Point& Point::operator=(Point another) {
	x = another.x;
	y = another.y;
	return *this;
}
//меньше та точка, которая левее и ниже
bool Point::operator<(Point another) const {
	return ((another.x - x) > 10e-7) || (comp(another.x, x) && ((another.y - y) > 10e-7));
}
//проверка совпадения
bool Point::operator==(Point another) const {
	return comp(x, another.x) && comp(y, another.y);
}
//проверка несовпадения
bool Point::operator!=(Point another) const {
	return !(*this == another);
}
//поворот
void Point::rotate(Point center, double angle) {
	double s = std::sin(M_PI * angle / 180);
	double c = std::cos(M_PI * angle / 180);
	double x_ = x;
	double y_ = y;
	x = center.x + c * (x_ - center.x) - s * (y_ - center.y);
	y = center.y + s * (x_ - center.x) + c * (y_ - center.y);
}
//центральная симметрия
void Point::reflex(Point center) {
	x = 2 * center.x - x;
	y = 2 * center.y - y;
}
//осевая симметрия
void Point::reflex(Line axis) {
	Point projection = Projection(*this, axis);
	reflex(projection);
}
//гомотетия
void Point::scale(Point center, double coefficient) {
	*this = center + Vector(center, *this) * coefficient;
}
//параллельный перенос (с созданием копии)
Point operator+(Point start, Vector transfer) {
	start += transfer;
	return start;
}
//середина отрезка
Point Middle(Point first, Point second) {
	return Point((Vector(first) + Vector(second)) / 2);
}
//точка, делящая отрезок в заданном отношении
Point Divide(Point first, Point second, double k) {
	k /= k + 1;
	Point divider = first + Vector(first, second) * k;
	return divider;
}
//расстояние между двумя точками
double Distance(Point first, Point second) {
	Vector v = Vector(first, second);
	return v.Len();
}
//проекция точки на прямую
Point Projection(Point point, Line line) {
	Line perpendikular = Per(point, line);
	Point projection = Intersect(perpendikular, line);
	return projection;
}
//расстояние от точки до прямой
double Distance(Point point, Line line) {
	return Distance(point, Projection(point, line));
}
//положение точки относительно прямой
bool Highground(Point point, Line line) {
	return (line.drct.x * (line.drct.x * (point.y - line.rds.x) - line.drct.y * (point.x - line.rds.x))) >= 0;
}








Line::Line(void) = default;
//прямая по точке и направляющему вектору
Line::Line(Point point, Vector direction) : rds(point), drct(direction) {}

//прямая по радиус-вектору и направляющему вектору
Line::Line(Vector radius, Vector direction) : rds(radius), drct(direction) {}

//прямая по точке и параллельной прямой
Line::Line(Point point, Line parallel) : rds(point), drct(parallel.drct) {}

//прямая по двум точкам
Line::Line(Point first, Point second) {
	rds = first;
	drct = Vector(first, second);
}
//прямая по угловому коэффициенту и сдвигу
Line::Line(double k, double b) {
	rds = Vector(0, b);
	drct = Vector(1, k);
}
//прямая по точке и угловому коэффициенту
Line::Line(Point point, double k) {
	rds = point;
	drct = Vector(1, k);
}
//проверка совпадения
bool Line::operator==(Line another) const {
	return (drct || another.drct) && (drct || Vector(rds - another.rds));
}
//проверка несовпадения
bool Line::operator!=(Line another) const {
	return !(*this == another);
}
//проверка параллельности
bool Line::operator||(Line another) const {
	return (drct || another.drct);
}
//поворот
void Line::rotate(Point center, double angle) {
	Point begin = Point(rds);
	begin.rotate(center, angle);
	rds = Vector(begin);
	drct.rotate(center, angle);
}
//центральная симметрия
void Line::reflex(Point center) {
	Point begin = Point(rds);
	begin.reflex(center);
	rds = Vector(begin);
}
//осевая симметрия
void Line::reflex(Line axis) {
	Point begin = Point(rds);
	begin.reflex(axis);
	rds = Vector(begin);
	drct.reflex(axis);
}
//гомотетия
void Line::scale(Point center, double coefficient) {
	Point begin = Point(rds);
	begin.scale(center, coefficient);
	rds = Vector(begin);
}
//пересечение прямых
Point Intersect(Line first, Line second) {
	if (first || second)
		return Point(INFINITY, INFINITY);
	else {
		double t = (second.drct ^ (first.rds - second.rds)) / (first.drct ^ second.drct);
		Point point = Point(first.rds + first.drct * t);
		return point;
	}
}
//перпендикуляр из точки на прямую
Line Per(Point point, Line line) {
	Vector direction = Vector(-line.drct.y, line.drct.x);
	return Line(point, direction);
}
//серединный перпендикуляр
Line MidPer(Point first, Point second) {
	return Per(Middle(first, second), Line(first, second));
}
//биссектриса
Line Bisector(Point A, Point O, Point B) {
	double AO = Distance(A, O);
	double BO = Distance(B, O);
	Point L = Divide(A, B, AO / BO);
	return Line(O, L);
}








Ellipse::Ellipse(void) = default;

Ellipse::Ellipse(Point F1_, Point F2_, double focal_distance_) :
	F1(F1_), F2(F2_), focal_distance(focal_distance_) {
	Standartization();
}
//фокусы
std::pair<Point, Point> Ellipse::focuses() const {
	return std::pair<Point, Point> {F1, F2};
}
//директрисы
std::pair<Line, Line> Ellipse::directrices() const {
	double a = focal_distance / 2;
	double c = Distance(F1, F2) / 2;
	double k = (a * a) / (c * c);
	Line l = Line(F1, F2);
	Vector v1 = Vector(center(), F1) * k;
	Vector v2 = Vector(center(), F2) * k;
	Point p1 = center() + v1;
	Point p2 = center() + v2;
	Line d1 = Per(p1, l);
	Line d2 = Per(p2, l);
	return std::pair<Line, Line>(d1, d2);
}
//эксцентриситет
double Ellipse::eccentricity() const {
	double c = Distance(F1, F2) / 2;
	double a = focal_distance / 2;
	return c / a;
}
//центр 
Point Ellipse::center() const {
	return Middle(F1, F2);
}
//F1 - левый, F2 - правый
void Ellipse::Standartization() {
	Point copy = F1;
	if (F2 < F1) {
		F1 = F2;
		F2 = copy;
	}
}
//периметр
double Ellipse::perimeter() const {
	double a = focal_distance / 2;
	double c = Distance(F1, F2) / 2;
	double b = std::sqrt(a * a - c * c);
	double ans = M_PI * (3 * (a + b) - std::sqrt((3 * a + b) * (a + 3 * b)));
	return ans;
}
//площадь
double Ellipse::area() const {
	double a = focal_distance / 2;
	double c = Distance(F1, F2) / 2;
	double b = std::sqrt(a * a - c * c);
	return M_PI * a * b;
}
//равенство (по множествам точек)
bool Ellipse::operator==(const Shape & another) const {
	bool is_Ellipse = dynamic_cast<const Ellipse*>(&another);
	if (!is_Ellipse)
		return false;
	const Ellipse& anotherEllipse = dynamic_cast<const Ellipse&>(another);
	return (((F1 == anotherEllipse.F1) && (F2 == anotherEllipse.F2)) || ((F2 == anotherEllipse.F1) && (F1 == anotherEllipse.F2))) &&
		comp(focal_distance, anotherEllipse.focal_distance);
}
//неравенство (по множествам точек)
bool Ellipse::operator!=(const Shape & another) const {
	return !(*this == another);
}
//равенство
bool Ellipse::isCongruentTo(const Shape & another) const {
	return Ellipse::isSimilarTo(another) && comp(Ellipse::area(), another.area());
}
//подобие
bool Ellipse::isSimilarTo(const Shape & another) const {
	bool is_Ellipse = dynamic_cast<const Ellipse*>(&another);
	if (!is_Ellipse)
		return false;
	const Ellipse& anotherEllipse = dynamic_cast<const Ellipse&>(another);
	return comp(Distance(F1, F2) / focal_distance,
		Distance(anotherEllipse.F1, anotherEllipse.F2) / anotherEllipse.focal_distance);
}
//содержит ли точку
bool Ellipse::containsPoint(Point point) const {
	return (Distance(F1, point) + Distance(F2, point) - focal_distance) < 1e-7;
}
//поворот
void Ellipse::rotate(Point center, double angle) {
	F1.rotate(center, angle);
	F2.rotate(center, angle);
	Standartization();
}
//центральная симметрия
void Ellipse::reflex(Point center) {
	F1.reflex(center);
	F2.reflex(center);
	Standartization();
}
//осевая симметрия
void Ellipse::reflex(Line axis) {
	F1.reflex(axis);
	F2.reflex(axis);
	Standartization();
}
//гомотетия
void Ellipse::scale(Point center, double coefficient) {
	F1.scale(center, coefficient);
	F2.scale(center, coefficient);
	focal_distance *= abs(coefficient);
	Standartization();
}








Circle::Circle(void) = default;
//окружность по центру и радиусу
Circle::Circle(Point center_, double radius_) :
	Ellipse(center_, center_, radius_ * 2), O(center_), R(radius_) {}
//центр
double Circle::radius() const {
	return R;
}
//периметр
double Circle::perimeter() const {
	return M_PI * R * 2;
}
//площадь
double Circle::area() const {
	return M_PI * R * R;
}
//поворот
void Circle::rotate(Point center, double angle) {
	Ellipse::rotate(center, angle);
	O = F1;
}
//центральная симметрия
void Circle::reflex(Point center) {
	Ellipse::reflex(center);
	O = F1;
}
//осевая симметрия
void Circle::reflex(Line axis) {
	Ellipse::reflex(axis);
	O = F1;
}
//гомотетия
void Circle::scale(Point center, double coefficient) {
	Ellipse::scale(center, coefficient);
	O = F1;
	R = focal_distance / 2;
}








Polygon::Polygon(void) = default;

Polygon::Polygon(const std::vector<Point> points_) : points(points_) {}
Polygon::Polygon(std::initializer_list<Point> list) {
	points.reserve(list.size());
	for (const Point& point : list)
		points.push_back(point);
}
//количество вершин
int Polygon::verticesCount() const {
	return points.size();
}
//вектор вершин
const std::vector<Point>& Polygon::getVertices() const {
	const std::vector<Point>& ans = points;
	return ans;
}
//отношение соседних сторон
double Polygon::ratio_of_side(size_t i) const {
	i = (i + points.size()) % points.size();
	size_t prev_i = (i - 1 + points.size()) % points.size();
	size_t next_i = (i + 1) % points.size();
	return (Distance(points[prev_i], points[i]) / Distance(points[next_i], points[i]));
}
//синус угла
double Polygon::angle(size_t i) const {
	i = (i + points.size()) % points.size();
	size_t prev_i = (i - 1 + points.size()) % points.size();
	size_t next_i = (i + 1) % points.size();
	Vector side1 = Vector(points[prev_i], points[i]);
	Vector side2 = Vector(points[next_i], points[i]);
	return ((side1 ^ side2) / (side1.Len() * side2.Len()));
}
//выпуклость
bool Polygon::isConvex() const {
	bool is_convex;
	for (size_t i = 0; i != points.size(); i++) {
		size_t j = (i + 1) % points.size();
		size_t k = (i - 1 + points.size()) % points.size();
		Vector v1 = Vector(points[i], points[k]);
		Vector v2 = Vector(points[i], points[j]);
		if (i == 0)
			is_convex = ((v1 ^ v2) < 0);
		if (is_convex != ((v1 ^ v2) < 0))
			return false;
	}
	return true;
}
//периметр
double Polygon::perimeter() const {
	double ans = 0;
	for (size_t i = 0; i != points.size(); i++) {
		size_t j = (i + 1) % points.size();
		ans += Distance(points[i], points[j]);
	}
	return ans;
}
//площадь
double Polygon::area() const {
	double ans = 0;
	for (size_t i = 1; i != points.size() - 1; i++) {
		Vector side = Vector(points[i], points[i + 1]);
		Vector diagonal = Vector(points[0], points[i]);
		ans += (side ^ diagonal) / 2;
	}
	return ans;
}
//равенство (по множествам точек)
bool Polygon::operator==(const Shape & another) const {
	bool is_Polygon = dynamic_cast<const Polygon*>(&another);
	if (!is_Polygon)
		return false;
	const Polygon& anotherPolygon = dynamic_cast<const Polygon&>(another);
	if (points.size() != anotherPolygon.points.size())
		return false;
	for (size_t i = 0; i != points.size(); i++) {
		for (size_t j = 0; j != points.size(); j++) {

			if (points[j] != anotherPolygon.points[(i + j) % points.size()])
				break;
			if (j != points.size() - 1)
				continue;
			return true;
		}
		for (size_t j = 0; j != points.size(); j++) {

			if (points[j] != anotherPolygon.points[(i - j + points.size()) % points.size()])
				break;
			if (j != points.size() - 1)
				continue;
			return true;
		}
	}
	return false;
}
//неравенство (по множествам точек)
bool Polygon::operator!=(const Shape & another) const {
	return !(*this == another);
}
//равенство
bool Polygon::isCongruentTo(const Shape & another) const {
	return Polygon::isSimilarTo(another) && comp(Polygon::area(), another.area());
}
//подобие
bool Polygon::isSimilarTo(const Shape & another) const {
	bool is_Polygon = dynamic_cast<const Polygon*>(&another);
	if (!is_Polygon)
		return false;
	const Polygon& anotherPolygon = dynamic_cast<const Polygon&>(another);
	if (points.size() != anotherPolygon.points.size())
		return false;
	for (size_t i = 0; i != points.size(); i++) {
		int x = ((Polygon::angle(0) * anotherPolygon.angle(i) < 0) ? -1 : 1);
		for (size_t j = 0; j != points.size(); j++) {
			bool is_similar = (comp(Polygon::ratio_of_side(j), anotherPolygon.ratio_of_side(i + j)) &&
				comp(x * Polygon::angle(j), anotherPolygon.angle(i + j)));
			if (!is_similar)
				break;
			if (j != (points.size() - 1))
				continue;
			return true;
		}
		for (size_t j = 0; j != points.size(); j++) {
			bool is_similar = (comp(Polygon::ratio_of_side(j), 1 / anotherPolygon.ratio_of_side(i - j)) &&
				comp(x * Polygon::angle(j), anotherPolygon.angle(i - j)));
			if (!is_similar)
				break;
			if (j != (points.size() - 1))
				continue;
			return true;
		}
	}
	return false;
}
//содержит ли точку
bool Polygon::containsPoint(Point point) const {
	Line line = Line(point, Vector(1, 0));
	int counter = 0;
	for (size_t i = 0; i != points.size(); i++) {
		size_t j = (i + 1) % points.size();
		Point Inter = Intersect(Line(points[i], points[j]), line);
		if ((Inter.x < point.x) && (((Inter.x - points[i].x) * (Inter.x - points[j].x)) < 1e-7)) {
			counter++;
		}
	}
	return (counter % 2);
}
//поворот
void Polygon::rotate(Point center, double angle) {
	for (size_t i = 0; i != points.size(); i++)
		points[i].rotate(center, angle);
}
//центральная симметрия
void Polygon::reflex(Point center) {
	for (size_t i = 0; i != points.size(); i++)
		points[i].reflex(center);
}
//осевая симметрия
void Polygon::reflex(Line axis) {
	for (size_t i = 0; i != points.size(); i++)
		points[i].reflex(axis);
}
//гомотетия
void Polygon::scale(Point center, double coefficient) {
	for (size_t i = 0; i != points.size(); i++)
		points[i].scale(center, coefficient);
}








Rectangle::Rectangle(void) = default;

Rectangle::Rectangle(Point first_, Point second_, double side_ratio_) :
	first(first_), second(second_), side_ratio(side_ratio_) {
	Point X = first;
	Point Y = second;
	X.rotate(center(), 90);
	Y.rotate(center(), 90);
	if (side_ratio > 1)
		side_ratio = 1 / side_ratio;
	if (side_ratio == 1) {
		points = { first, Y, second, X };
		Standartization();
	}
	else {
		Point L1 = Divide(first, second, side_ratio);
		Point L2 = Divide(first, second, -side_ratio);
		Point third = Intersect(Line(X, L1), Line(Y, L2));
		Point fourth = third;
		fourth.reflex(center());
		points = { first, third, second, fourth };
		Standartization();
	}
}
//центр
Point Rectangle::center() const {
	return Middle(first, second);
}
//диагонали
std::pair<Line, Line> Rectangle::diagonals() const {
	Line diagonal1 = Line(points[0], points[2]);
	Line diagonal2 = Line(points[1], points[3]);
	return std::pair<Line, Line> {diagonal1, diagonal2};
}
//переобозначение вершин
void Rectangle::Standartization() {
	first = points[0];
	second = points[2];
}
//поворот
void Rectangle::rotate(Point center, double angle) {
	Polygon::rotate(center, angle);
	Standartization();
}
//центральная симметрия
void Rectangle::reflex(Point center) {
	Polygon::reflex(center);
	Standartization();
}
//осевая симметрия
void Rectangle::reflex(Line axis) {
	Polygon::reflex(axis);
	Standartization();
}
//гомотетия
void Rectangle::scale(Point center, double coefficient) {
	Polygon::scale(center, coefficient);
	Standartization();
}








Square::Square(void) = default;

Square::Square(Point first_, Point second_) : Rectangle(first_, second_, 1) {}
//описанная окружность
Circle Square::circumscribedCircle() const {
	Point center = Middle(first, second);
	double radius = Distance(first, second) / 2;
	return Circle(center, radius);
}
//вписанная окружность
Circle Square::inscribedCircle() const {
	Point center = Middle(first, second);
	double radius = Distance(first, second) / (2 * std::sqrt(2));
	return Circle(center, radius);
}








Triangle::Triangle(void) = default;

Triangle::Triangle(Point A_, Point B_, Point C_) : Polygon({ A_, B_, C_ }), A(A_), B(B_), C(C_) {}
//описанная окружность
Circle Triangle::circumscribedCircle() const {
	Point O = circumcenter();
	double R = Distance(A, O);
	return Circle(O, R);
}
//вписанная окружность
Circle Triangle::inscribedCircle() const {
	Point I = incenter();
	double r = Distance(I, Line(B, C));
	return Circle(I, r);
}
//центр описанной окружности
Point Triangle::circumcenter() const {
	Point O = Intersect(MidPer(A, B), MidPer(A, C));
	return O;
}
//инцентр
Point Triangle::incenter() const {
	Point I = Intersect(Bisector(A, B, C), Bisector(A, C, B));
	return I;
}
//центроид
Point Triangle::centroid() const {
	Point G = Point((A.x + B.x + C.x) / 3, (A.y + B.y + C.y) / 3);
	return G;
}
//ортоцентр
Point Triangle::orthocenter() const {
	Point O = circumcenter();
	Point M = Middle(B, C);
	Point H = A + Vector(O, M) * 2;
	return H;
}
//прямая Эйлера
Line Triangle::EulerLine() const {
	Point O = circumcenter();
	Point H = orthocenter();
	Line OH = Line(O, H);
	return OH;
}
//окружность Эйлера/9-ти точек
Circle Triangle::ninePointsCircle() const {
	Point O = circumcenter();
	Point H = orthocenter();
	Point X = Middle(O, H);
	double R = circumscribedCircle().R;
	return Circle(X, R / 2);
}
//переобозначение вершин
void Triangle::Standartization() {
	A = points[0];
	B = points[1];
	C = points[2];
}
//поворот
void Triangle::rotate(Point center, double angle) {
	Polygon::rotate(center, angle);
	Standartization();
}
//центральная симметрия
void Triangle::reflex(Point center) {
	Polygon::reflex(center);
	Standartization();
}
//осевая симметрия
void Triangle::reflex(Line axis) {
	Polygon::reflex(axis);
	Standartization();
}
//гомотетия
void Triangle::scale(Point center, double coefficient) {
	Polygon::scale(center, coefficient);
	Standartization();
}
