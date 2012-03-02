#include <string.h>
#include "support.h"


// -- class Point2D ---
Point2D::Point2D () : x_ (0),
		      y_ (0) {
}

Point2D::Point2D (float x, float y) : x_ (x),
			 	      y_ (y) {
}

Point2D::~Point2D () {

}

Point2D& Point2D::operator= (const Point2D& point) {
    if (this == &point)
        return *this;
    x_ = point.x_;
    y_ = point.y_;

    return *this;
}
