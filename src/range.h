
#ifndef _range_H_
#define _range_H_

template<class T>
class general_range {
 public:
   class iterator {
      friend class general_range;
    public:
      T operator *() const { return i_; }
      iterator& operator ++() { ++i_; return *this; }

      bool operator ==(const iterator &other) const { return i_ == other.i_; }
      bool operator !=(const iterator &other) const { return i_ != other.i_; }

    protected:
      iterator(T start) : i_ (start) { }

    private:
      T i_;
   };

   iterator begin() const { return begin_; }
   iterator end() const { return end_; }
   general_range(T  begin, T end) : begin_(begin), end_(end) 
   { if (end < begin) end = begin; }
private:
   iterator begin_;
   iterator end_;
};

inline
general_range<long> range(long n) { return general_range<long>(0, n); }

inline
general_range<long> range(long m, long n) { return general_range<long>(m, n); }

#endif

