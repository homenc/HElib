/* Copyright (C) 2012,2013 IBM Corp.
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
#ifndef CLONED_PTR_H
#define CLONED_PTR_H
/**
* @file cloned_ptr.h
* @brief Implemenation of smart pointers with "deep cloning" semantics
*
* Based (loosely) on code from
*
*  http://github.com/yonat/smart_ptr/blob/master/cloned_ptr.h
**/
/*
* This seems like a very useful template class for smart
* pointers with "deep cloning" semantics and class heirarchies.
* 
* To work, a class needs a "clone" method that creates a new copy.
* 
* NOTES:
* 
* this did not compile...had to change the definition of the class clone<X> 
* 
* added boolean operators
* 
* added assign operation set_ptr(X* p) removed get(), introduced public
* method get_ptr()... these names chosen to make it easier to hunt down
* "loopholes"
* 
* added noclone<X> template
* 
* Made "const"-ness apply to both the pointer and the object pointed to.
* For this kind of "deep copy" semantics, this makes the most sense.
* 
* Made a corresponding class copied_ptr, with "shallow" cloning...
* To bad this can't really be done with templates... C++11 defines
* the concept of "template aliases", which would do the trick...
* unfortunately, gcc does not yet implement this...
* 
* Changed the name "clone" to avoid conflicts...
*****/

/* For ANSI-challenged compilers, you may want to
 * #define NO_MEMBER_TEMPLATES or explicit */

/**
* @class deep_clone
* @brief Deep copy: initialize with clone
* @tparam X The class to which this points
**/
template <class X> class deep_clone
{
public:
    static X* apply(const X* x) {return x->clone();}
};

/**
* @class shallow_clone
* @brief Shallow copy: initialize with copy constructor
* @tparam X The class to which this points
**/
template <class X> class shallow_clone
{
public:
    static X* apply(const X* x) {return new X(*x); }
};


#ifndef NO_MEMBER_TEMPLATES

#define CLONED_PTR_TEMPLATE_MEMBERS(CLONED_PTR_TYPE) \
 \
    template <class Y> CLONED_PTR_TYPE(const CLONED_PTR_TYPE<Y>& r) \
        {copy(r.ptr);} \
    template <class Y> CLONED_PTR_TYPE& operator=(const CLONED_PTR_TYPE<Y>& r) \
    { \
        if (this != &r) { \
            delete ptr; \
            copy(r.ptr); \
        } \
        return *this; \
    } \

#else

#define CLONED_PTR_TEMPLATE_MEMBERS(CLONED_PTR_TYPE) 

#endif



#define CLONED_PTR_DECLARE(CLONED_PTR_TYPE,CLONED_PTR_INIT) \
 \
template <class X, class Cloner = CLONED_PTR_INIT<X> > class CLONED_PTR_TYPE \
{ \
public: \
    typedef X element_type; \
 \
    explicit CLONED_PTR_TYPE(X* p = 0)  : ptr(p) {} \
    ~CLONED_PTR_TYPE() {delete ptr;} \
    CLONED_PTR_TYPE(const CLONED_PTR_TYPE& r) {copy(r.ptr);} \
 \
    CLONED_PTR_TYPE& operator=(const CLONED_PTR_TYPE& r) \
    { \
        if (this != &r) { \
            delete ptr; \
            copy(r.ptr); \
        } \
        return *this; \
    } \
 \
    void set_ptr(X* p) \
    { \
        if (ptr != p) { \
            delete ptr; \
            ptr = p; \
        } \
    } \
 \
    CLONED_PTR_TEMPLATE_MEMBERS(CLONED_PTR_TYPE) \
 \
    const X& operator*()  const  {return *ptr;} \
    X& operator*() {return *ptr;} \
 \
    const X* operator->() const {return ptr;} \
    X* operator->() {return ptr;} \
 \
    bool null() const { return ptr == 0; } \
 \
    const X* get_ptr() const { return ptr; } \
    X* get_ptr() { return ptr; } \
 \
    void swap(CLONED_PTR_TYPE& r) \
    { \
      X *tmp; tmp = r.ptr; r.ptr = ptr; ptr = tmp;  \
    } \
private: \
    X* ptr; \
 \
    void copy(X* p)  {ptr = (p ? Cloner::apply(p) : 0);} \
}; \
 \



// declare the template class cloned_ptr<X>
CLONED_PTR_DECLARE(cloned_ptr, deep_clone)

// declare the template class coppied_ptr<X>
CLONED_PTR_DECLARE(copied_ptr, shallow_clone)

// template <class X> using copied_ptr = cloned_ptr<X, noclone<X> >;

template<class X, class Cloner> 
void swap(cloned_ptr<X,Cloner>& x, cloned_ptr<X,Cloner>& y) { x.swap(y); }

template<class X, class Cloner> 
void swap(copied_ptr<X,Cloner>& x, copied_ptr<X,Cloner>& y) { x.swap(y); }



#endif // CLONED_PTR_H

