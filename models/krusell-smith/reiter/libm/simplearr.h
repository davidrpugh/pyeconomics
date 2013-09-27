#ifndef FM_SISIMPLEARR_H_
#define FM_SISIMPLEARR_H_ 1

#include <assert.h>

template<class T>
class SimpleArray {
  T *x;
  size_t nr,nc;
  int owns_data;
 public:
  SimpleArray(){nr=0;nc=0;owns_data=0;}
  SimpleArray(size_t n1,size_t n2)
    {
      owns_data=0;
      resize(n1,n2);

    }
  SimpleArray(size_t n1)
    {
      owns_data=0;
      resize(n1,1);
    }
  SimpleArray(T *x0, size_t n1)
    {
      owns_data=0;
      x = x0;
      nr = n1;
      nc = 1;
    }
  SimpleArray(T *x0, size_t n1, size_t n2)
    {
      owns_data=0;
      x = x0;
      nr = n1;
      nc = n2;
    }
  ~SimpleArray() {if(owns_data) delete[] x;owns_data=0;}
  
  void resize(size_t n1,size_t n2) 
  {
    if(owns_data) delete[] x;
    nr=n1;nc=n2;
    if(nr*nc>0){
      x = new T[nr*nc];
      owns_data=1;
    } else 
      owns_data=0;
  }
  void resize(size_t n1) {resize(n1,1);}
  
  SimpleArray(const SimpleArray &source){owns_data=0;*this = source;}
  SimpleArray &operator=(const SimpleArray<T>&source)
    {
      if( (void*)this == (void*) &source)
	return *this;
      resize(source.nr,source.nc);
      size_t i,n=size();
      for(i=0;i<n;i++)
	x[i] = source.x[i];
      return *this;
    };
  int Owns() const {return owns_data;}
  void SetNew() {owns_data=0;}
  T &operator[](size_t i){
    myassert(i<size(),"subscript error");
    return x[i];
  }
  T &operator[](size_t i) const {
    myassert(i<size(),"subscript error");
    return x[i];
  }
  T &operator()(size_t i){
    i--;
    myassert(i<size(),"subscript error");
    return (*this)[i];
  }
  T &operator()(size_t i) const {
    i--;
    myassert(i<size(),"subscript error");
    return (*this)[i];
  }

  T &operator()(size_t i,size_t j){
    i--;
    j--;
    myassert(i<nr && j<nc,"subscript error");
    return x[j*nr + i];
  }
  T &operator()(size_t i,size_t j) const {
    i--;
    j--;
    myassert(i<nr && j<nc,"subscript error");
    return x[j*nr + i];
  }

  size_t size() const {return nr*nc;}
  size_t size(size_t i) const {myassert(i==1||i==2,"subscript error");return i==1 ? nr : nc;}
  T *begin() {return x;};
  const T *begin() const {return x;};
  T *end() {return x+size();};
  const T *end() const {return x+size();};
};


#endif
