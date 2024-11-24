#ifndef VEC_HH_INCLUDED
#define VEC_HH_INCLUDED

#include <cmath>
#include <iostream>
#include <sstream>

/**
 * The Vec class is a straightforward 3-d vector, which also works in pyroot.
 */
struct Vec
{
  double x,y,z;

  /**
   * Constructor.
   *
   * \param  x_            x position
   * \param  y_            y position
   * \param  z_            z position
   */
  Vec(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
  
  /**
   * Default constructor.
   */
  Vec():x(0),y(0),z(0) {}
  
  /**
   * Get dot product.
   *
   * \param  v             vector
   * \return               dot product
   */
  double dot(const Vec& v) const { return v.x*x + v.y*y+ v.z*z;}
  
  /**
   * Get cross product.
   *
   * \param  r             vector
   * \return               cross product
   */
  Vec cross(const Vec r)   const { return Vec ( y*r.z-z*r.y, z*r.x-x*r.z, x*r.y-y*r.x);}

  /**
   * Add vector.
   *
   * \param  v             vector
   * \return               this vector
   */
  Vec& operator+=(const Vec& v)  { x+=v.x; y+=v.y; z+=v.z; return *this;}

  /**
   * Subtract vector.
   *
   * \param  v             vector
   * \return               this vector
   */
  Vec& operator-=(const Vec& v)  { x-=v.x; y-=v.y; z-=v.z; return *this;}

  /**
   * Multiply vector.
   *
   * \param  d             factor
   * \return               this vector
   */
  Vec& operator*=(double d)      { x*=d; y*=d; z*=d; return *this;}

  /**
   * Divide vector.
   *
   * \param  d             factor
   * \return               this vector
   */
  Vec& operator/=(double d)      { return operator*=( 1.0 / d ); }

  /**
   * Negate vector.
   *
   * \return               vector
   */
  Vec operator-() const          { return Vec(-x,-y,-z); }
  
  /**
   * Check equality with given vector.
   *
   * \param  v             vector
   * \return               true if (x,y,z) positions of two vectors are equal; else false
   */
  bool operator==( const Vec& v ) const { return x==v.x && y==v.y && z==v.z ; }

  /**
   * Check in-equality with given vector.
   *
   * \param  v             vector
   * \return               true if one of (x,y,z) positions of two vectors are not equal; else false
   */
  bool operator!=( const Vec& v ) const { return x!=v.x || y!=v.y || z!=v.z ; }

  /**
   * Set vector.
   *
   * \param  xx            x position
   * \param  yy            y position
   * \param  zz            z position
   * \return               this vector
   */
  Vec& set(double xx, double yy, double zz) { x=xx; y=yy; z=zz; return *this;}

  /**
   * Set vector according given zenith and azimuth angles.
   *
   * \param  theta         zenith  angle [rad]
   * \param  phi           azimuth angle [rad]
   * \return               this vector
   */
  Vec& set_angles(double theta, double phi)
  {
    x = sin ( theta ) * cos( phi );
    y = sin ( theta ) * sin( phi );
    z = cos ( theta );
    return *this;
  }

  Vec& set_spherical(double r, double theta, double phi)
  {
    x = r * sin ( theta ) * cos( phi );
    y = r * sin ( theta ) * sin( phi );
    z = r * cos ( theta );
    return *this;
  }

  /**
   * Get azimuth angle.
   *
   * \return               angle [rad]
   */
  double phi()   const { return atan2( y,x ); }

  /**
   * Get zenith angle.
   *
   * \return               angle [rad]
   */
  double theta() const { return acos(z / len()); }

  /**
   * Get length.
   *
   * \return               length
   */
  double len()   const { double l = dot(*this); return (l > 0)? sqrt(l) : 0; }

  /**
   * Get length of (x,y) component.
   *
   * \return               length
   */
  double lenxy() const { const double r2 = x*x + y*y; return (r2>0) ? sqrt(r2) :0; }

  /**
   * Normalise this vector.
   *
   * \return               this vector
   */
  Vec& normalize() { return operator/=( len() ); }

  /**
   * Print vector.
   *
   * \param  out           output stream
   */
  void print( std::ostream& out = std::cout ) const
  {
  	out << "Vec:" << x << " " << y << " " << z;
  }

  /**
   * Get string representation of this vector
   *
   * \return               string
   */
  const char* __repr__() const 
  {
    static std::string buffer;

    std::ostringstream s; 
  
    print(s);

    buffer = s.str();

    return buffer.c_str();
  }

  /**
   * Add vector.
   *
   * \param  v             vector
   * \return               vector
   */
  Vec __add__(const Vec& v) const { Vec r=*this; return r+=v; }

  /**
   * Subtract vector.
   *
   * \param  v             vector
   * \return               vector
   */
  Vec __sub__(const Vec& v) const { Vec r=*this; return r-=v; }

  /**
   * Multiply vector.
   *
   * \param  d             factor
   * \return               vector
   */
  Vec __mul__(double d )    const { Vec r=*this; return r*=d; }

  /**
   * Multiply vector.
   *
   * \param  d             factor
   * \return               vector
   */
  Vec __rmul__(double d )   const { return __mul__(d);        } 

  /**
   * Divide vector.
   *
   * \param  d             factor
   * \return               vector
   */
  Vec __div__(double d )    const { Vec r=*this; return r/=d; }

  /**
   * Rotate around z-axis with given angle.
   *
   * \param  ang           angle [rad]
   * \return               this vector
   */
  Vec& rotate_z(double ang)
  {
    const Vec o = *this; 
    x = o.x *cos(ang) - o.y * sin(ang);
    y = o.x *sin(ang) + o.y * cos(ang);
    z = o.z;
    return *this;
  }

  /**
   * Rotate around x-axis with given angle.
   *
   * \param  ang           angle [rad]
   * \return               this vector
   */
  Vec& rotate_x(double ang)
  {
    const Vec o = *this;
    x = o.x;
    y = o.y *cos(ang) + o.z * -sin(ang);
    z = o.y *sin(ang) + o.z * cos(ang);
    return *this;
  }

  /**
   * Rotate around y-axis with given angle.
   *
   * \param  ang           angle [rad]
   * \return               this vector
   */
  Vec& rotate_y(double ang)
  {
    const Vec o = *this;
    x = o.x *cos(ang) + o.z * sin(ang);
    y = o.y;
    z = -o.x *sin(ang) + o.z * cos(ang);
    return *this;
  }
};

/**
 * Write vector to output stream.
 *
 * \param  out           output stream
 * \param  v             vector
 * \return               output stream
 */
inline std::ostream& operator<<( std::ostream& out , const Vec& v )
{
  out << v.x << " " << v.y << " " << v.z << " ";
  return out;
}

/**
 * Read vector from input stream.
 *
 * \param  in            input stream
 * \param  v             vector
 * \return               input stream
 */
inline std::istream& operator>>(std::istream& in, Vec& v) 
{ 
  in >> v.x >> v.y >> v.z ; return in;
}

/**
 * Get cosine of space angle between two vectors.
 *
 * \param  a             first  vector
 * \param  b             second vector
 * \return               cosine
 */
inline double cos_angle_between( const Vec& a, const Vec& b)
{
  const double n = a.len() * b.len();
  return a.dot(b) / n;
}

/**
 * Get space angle between two vectors.
 *
 * \param  a             first  vector
 * \param  b             second vector
 * \return               angle [rad]
 */
inline double angle_between( const Vec& a, const Vec& b )
{
  double c = cos_angle_between( a, b );
  if ( c < -1 ) return M_PI;
  if ( c > 1  ) return 0;
  return acos( c );
}

/**
 * Add two vectors.
 *
 * \param  a             first  vector
 * \param  b             second vector
 * \return               vector
 */
inline Vec operator+(const Vec& a, const Vec& b) { Vec r(a); return r+=b;}

/**
 * Subtract two vectors.
 *
 * \param  a             first  vector
 * \param  b             second vector
 * \return               vector
 */
inline Vec operator-(const Vec& a, const Vec& b) { Vec r(a); return r-=b;}

/**
 * Multiply vector.
 *
 * \param  a             factor
 * \param  v             vector
 * \return               vector
 */
inline Vec operator*(double a, const Vec& v)     { return Vec(a*v.x,a*v.y,a*v.z);}

/**
 * Multiply vector.
 *
 * \param  v             vector
 * \param  a             factor
 * \return               vector
 */
inline Vec operator*(const Vec& v, double a)     { return Vec(a*v.x,a*v.y,a*v.z);}

/**
 * Divide vector.
 *
 * \param  v             vector
 * \param  a             factor
 * \return               vector
 */
inline Vec operator/(const Vec& v, double a)     { return Vec(v.x/a,v.y/a,v.z/a);}

#endif
