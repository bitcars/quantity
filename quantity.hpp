//
// Created by Yiding Han on 5/8/21.
//

#ifndef UNIT_DATASIZE_HPP
#define UNIT_DATASIZE_HPP

#include <type_traits>
#include <ratio>
#include <limits>

namespace datasize
{

template <class Rep, class Scale = std::ratio<1> > class unit;

template <class Tp>
struct _is_unit : std::false_type {};

template<class Rep, class Scale>
struct _is_unit<unit<Rep, Scale> > : std::true_type  {};

template <class Rep, class Scale>
struct _is_unit<const unit<Rep, Scale> > : std::true_type  {};

template <class Rep, class Scale>
struct _is_unit<volatile unit<Rep, Scale> > : std::true_type  {};

template <class Rep, class Scale>
struct _is_unit<const volatile unit<Rep, Scale> > : std::true_type  {};

} // datasize

template <class Rep1, class Scale1, class Rep2, class Scale2>
struct std::common_type<datasize::unit<Rep1, Scale1>,
datasize::unit<Rep2, Scale2> >
{
typedef datasize::unit<typename std::common_type<Rep1, Rep2>::type,
        typename __ratio_gcd<Scale1, Scale2>::type> type;
};

namespace datasize
{

// unit_cast
template <class FromUnit, class ToUnit,
        class Scale = typename std::ratio_divide<typename FromUnit::scale, typename ToUnit::scale>::type,
        bool = Scale::num == 1,
        bool = Scale::den == 1>
struct _unit_cast;

template <class FromUnit, class ToUnit, class Scale>
struct _unit_cast<FromUnit, ToUnit, Scale, true, true>
{
    constexpr
    ToUnit operator()(const FromUnit& fu) const
    {
        return ToUnit(static_cast<typename ToUnit::rep>(fu.count()));
    }
};

template <class FromUnit, class ToUnit, class Scale>
struct _unit_cast<FromUnit, ToUnit, Scale, true, false>
{
    constexpr
    ToUnit operator()(const FromUnit& fu) const
    {
        typedef typename std::common_type<typename ToUnit::rep, typename FromUnit::rep, intmax_t>::type Ct;
        return ToUnit(static_cast<typename ToUnit::rep>(
                                   static_cast<Ct>(fu.count()) / static_cast<Ct>(Scale::den)));
    }
};

template <class FromUnit, class ToUnit, class Scale>
struct _unit_cast<FromUnit, ToUnit, Scale, false, true>
{
    constexpr
    ToUnit operator()(const FromUnit& fu) const
    {
        typedef typename std::common_type<typename ToUnit::rep, typename FromUnit::rep, intmax_t>::type Ct;
        return ToUnit(static_cast<typename ToUnit::rep>(
                                   static_cast<Ct>(fu.count()) * static_cast<Ct>(Scale::num)));
    }
};

template <class FromUnit, class ToUnit, class Scale>
struct _unit_cast<FromUnit, ToUnit, Scale, false, false>
{
    constexpr
    ToUnit operator()(const FromUnit& fu) const
    {
        typedef typename std::common_type<typename ToUnit::rep, typename FromUnit::rep, intmax_t>::type Ct;
        return ToUnit(static_cast<typename ToUnit::rep>(
                                   static_cast<Ct>(fu.count()) * static_cast<Ct>(Scale::num)
                                   / static_cast<Ct>(Scale::den)));
    }
};

template <class ToUnit, class Rep, class Scale>
inline constexpr
typename std::enable_if
    <
        _is_unit<ToUnit>::value,
        ToUnit
    >::type
unit_cast(const unit<Rep, Scale>& fu)
{
    return _unit_cast<unit<Rep, Scale>, ToUnit>()(fu);
}

template <class Rep>
struct treat_as_floating_point : std::is_floating_point<Rep> {};

template <class Rep>
struct unit_values
{
public:
    static constexpr Rep zero() noexcept {return Rep(0);}
    static constexpr Rep max()  noexcept {return std::numeric_limits<Rep>::max();}
    static constexpr Rep min()  noexcept {return std::numeric_limits<Rep>::lowest();}
};

template <class Rep, class Scale>
class unit
{
    static_assert(!_is_unit<Rep>::value, "A unit representation can not be a unit");
    static_assert(std::__is_ratio<Scale>::value, "Second template parameter of unit must be a std::ratio");
    static_assert(Scale::num > 0, "unit scale must be positive");

    template <class R1, class R2>
    struct _no_overflow
    {
    private:
        static const intmax_t gcd_n1_n2_ = std::__static_gcd<R1::num, R2::num>::value;
        static const intmax_t gcd_d1_d2_ = std::__static_gcd<R1::den, R2::den>::value;
        static const intmax_t n1_ = R1::num / gcd_n1_n2_;
        static const intmax_t d1_ = R1::den / gcd_d1_d2_;
        static const intmax_t n2_ = R2::num / gcd_n1_n2_;
        static const intmax_t d2_ = R2::den / gcd_d1_d2_;
        static const intmax_t max_ = -((intmax_t(1) << (sizeof(intmax_t) * CHAR_BIT - 1)) + 1);

        template <intmax_t Xp, intmax_t Yp, bool overflow>
        struct _mul    // overflow == false
        {
            static const intmax_t value = Xp * Yp;
        };

        template <intmax_t Xp, intmax_t Yp>
        struct _mul<Xp, Yp, true>
        {
            static const intmax_t value = 1;
        };

    public:
        static const bool value = (n1_ <= max_ / d2_) && (n2_ <= max_ / d1_);
        typedef std::ratio<_mul<n1_, d2_, !value>::value,
        _mul<n2_, d1_, !value>::value> type;
    };

public:
    typedef Rep rep;
    typedef typename Scale::type scale;
private:
    rep rep_;
public:

    constexpr unit() = default;

    template <class Rep2>
    constexpr
    explicit unit(const Rep2& r,
        typename std::enable_if
            <
                std::is_convertible<Rep2, rep>::value &&
                (treat_as_floating_point<rep>::value ||
                 !treat_as_floating_point<Rep2>::value)
            >::type* = nullptr)
        : rep_(r) {}

    // conversions

    template <class Rep2, class Scale2>
    constexpr
    unit(const unit<Rep2, Scale2>& u,
        typename std::enable_if
            <
                _no_overflow<Scale2, scale>::value && (
                    treat_as_floating_point<rep>::value ||
                    (_no_overflow<Scale2, scale>::type::den == 1 &&
                     !treat_as_floating_point<Rep2>::value))
            >::type* = nullptr)
        : rep_(datasize::unit_cast<unit>(u).count()) {}

    // observer

    constexpr rep count() const {return rep_;}

    // arithmetic

    constexpr typename std::common_type<unit>::type operator+() const {return typename std::common_type<unit>::type(*this);}
    constexpr typename std::common_type<unit>::type operator-() const {return typename std::common_type<unit>::type(-rep_);}
    unit& operator++()      {++rep_; return *this;}
    unit  operator++(int)   {return unit(rep_++);}
    unit& operator--()      {--rep_; return *this;}
    unit  operator--(int)   {return unit(rep_--);}

    unit& operator+=(const unit& u) {rep_ += u.count(); return *this;}
    unit& operator-=(const unit& u) {rep_ -= u.count(); return *this;}

    unit& operator*=(const rep& rhs) {rep_ *= rhs; return *this;}
    unit& operator/=(const rep& rhs) {rep_ /= rhs; return *this;}
    unit& operator%=(const rep& rhs) {rep_ %= rhs; return *this;}
    unit& operator%=(const unit& rhs) {rep_ %= rhs.count(); return *this;}

    // special values

    static constexpr unit zero() noexcept {return unit(unit_values<rep>::zero());}
    static constexpr unit min()  noexcept {return unit(unit_values<rep>::min());}
    static constexpr unit max()  noexcept {return unit(unit_values<rep>::max());}
};

typedef unit< long long,          std::ratio<1> > i_bytes;
typedef unit< long long,          std::ratio<1> > i_kilobytes;
typedef unit< long long,       std::ratio<1024> > i_megabytes;
typedef unit<      long,    std::ratio<1048576> > i_gigabytes;
typedef unit<      long, std::ratio<1073741824> > i_terabytes;

typedef unit< double,    std::ratio<1, 1024> > bytes;
typedef unit< double,          std::ratio<1> > kilobytes;
typedef unit< double,       std::ratio<1024> > megabytes;
typedef unit< double,    std::ratio<1048576> > gigabytes;
typedef unit< double, std::ratio<1073741824> > terabytes;


// Unit ==

template <class LhsUnit, class RhsUnit>
struct _unit_eq
{
    constexpr
    bool operator()(const LhsUnit& lhs, const RhsUnit& rhs) const
    {
        typedef typename std::common_type<LhsUnit, RhsUnit>::type Ct;
        return Ct(lhs).count() == Ct(rhs).count();
    }
};

template <class LhsUnit>
struct _unit_eq<LhsUnit, LhsUnit>
{
    constexpr
    bool operator()(const LhsUnit& lhs, const LhsUnit& rhs) const
    {return lhs.count() == rhs.count();}
};

template <class Rep1, class Scale1, class Rep2, class Scale2>
inline constexpr bool
operator==(const unit<Rep1, Scale1>& lhs, const unit<Rep2, Scale2>& rhs)
{
    return _unit_eq<unit<Rep1, Scale1>, unit<Rep2, Scale2> >()(lhs, rhs);
}

// Unit !=

template <class Rep1, class Scale1, class Rep2, class Scale2>
inline constexpr bool
operator!=(const unit<Rep1, Scale1>& lhs, const unit<Rep2, Scale2>& rhs)
{
    return !(lhs == rhs);
}

// Unit <

template <class LhsUnit, class RhsUnit>
struct _unit_lt
{
    constexpr
    bool operator()(const LhsUnit& lhs, const RhsUnit& rhs) const
    {
        typedef typename std::common_type<LhsUnit, RhsUnit>::type Ct;
        return Ct(lhs).count() < Ct(rhs).count();
    }
};

template <class LhsUnit>
struct _unit_lt<LhsUnit, LhsUnit>
{
    constexpr
    bool operator()(const LhsUnit& lhs, const LhsUnit& rhs) const
    {
        return lhs.count() < rhs.count();
    }
};

template <class Rep1, class Scale1, class Rep2, class Scale2>
inline constexpr bool
operator< (const unit<Rep1, Scale1>& lhs, const unit<Rep2, Scale2>& rhs)
{
    return _unit_lt<unit<Rep1, Scale1>, unit<Rep2, Scale2> >()(lhs, rhs);
}

// Unit >

template <class Rep1, class Scale1, class Rep2, class Scale2>
inline constexpr bool
operator> (const unit<Rep1, Scale1>& lhs, const unit<Rep2, Scale2>& rhs)
{
    return rhs < lhs;
}

// Unit <=

template <class Rep1, class Scale1, class Rep2, class Scale2>
inline constexpr bool
operator<=(const unit<Rep1, Scale1>& lhs, const unit<Rep2, Scale2>& rhs)
{
    return !(rhs < lhs);
}

// Unit >=

template <class Rep1, class Scale1, class Rep2, class Scale2>
inline constexpr bool
operator>=(const unit<Rep1, Scale1>& lhs, const unit<Rep2, Scale2>& rhs)
{
    return !(lhs < rhs);
}

// Unit +

template <class Rep1, class Scale1, class Rep2, class Scale2>
inline constexpr
typename std::common_type<unit<Rep1, Scale1>, unit<Rep2, Scale2> >::type
operator+(const unit<Rep1, Scale1>& lhs, const unit<Rep2, Scale2>& rhs)
{
    typedef typename std::common_type<unit<Rep1, Scale1>, unit<Rep2, Scale2> >::type Cu;
    return Cu(Cu(lhs).count() + Cu(rhs).count());
}

// Unit -

template <class Rep1, class Scale1, class Rep2, class Scale2>
inline constexpr
typename std::common_type<unit<Rep1, Scale1>, unit<Rep2, Scale2> >::type
operator-(const unit<Rep1, Scale1>& lhs, const unit<Rep2, Scale2>& rhs)
{
    typedef typename std::common_type<unit<Rep1, Scale1>, unit<Rep2, Scale2> >::type Cu;
    return Cu(Cu(lhs).count() - Cu(rhs).count());
}

// Unit *

template <class Rep1, class Scale, class Rep2>
inline constexpr
typename std::enable_if
    <
        std::is_convertible<Rep2, typename std::common_type<Rep1, Rep2>::type>::value,
        unit<typename std::common_type<Rep1, Rep2>::type, Scale>
    >::type
operator*(const unit<Rep1, Scale>& u, const Rep2& s)
{
    typedef typename std::common_type<Rep1, Rep2>::type Cr;
    typedef unit<Cr, Scale> Cu;
    return Cu(Cu(u).count() * static_cast<Cr>(s));
}

template <class Rep1, class Scale, class Rep2>
inline constexpr
typename std::enable_if
    <
        std::is_convertible<Rep1, typename std::common_type<Rep1, Rep2>::type>::value,
        unit<typename std::common_type<Rep1, Rep2>::type, Scale>
    >::type
operator*(const Rep1& s, const unit<Rep2, Scale>& u)
{
    return u * s;
}

// Unit /

template <class Rep1, class Scale, class Rep2>
inline constexpr
typename std::enable_if
    <
        !_is_unit<Rep2>::value &&
        std::is_convertible<Rep2, typename std::common_type<Rep1, Rep2>::type>::value,
        unit<typename std::common_type<Rep1, Rep2>::type, Scale>
    >::type
operator/(const unit<Rep1, Scale>& u, const Rep2& s)
{
    typedef typename std::common_type<Rep1, Rep2>::type Cr;
    typedef unit<Cr, Scale> Cu;
    return Cu(Cu(u).count() / static_cast<Cr>(s));
}

template <class Rep1, class Scale1, class Rep2, class Scale2>
inline constexpr typename std::common_type<Rep1, Rep2>::type
operator/(const unit<Rep1, Scale1>& lhs, const unit<Rep2, Scale2>& rhs)
{
    typedef typename std::common_type<unit<Rep1, Scale1>, unit<Rep2, Scale2> >::type Ct;
    return Ct(lhs).count() / Ct(rhs).count();
}

// Unit %

template <class Rep1, class Scale, class Rep2>
inline constexpr typename std::enable_if
    <
        !_is_unit<Rep2>::value &&
        std::is_convertible<Rep2, typename std::common_type<Rep1, Rep2>::type>::value,
        unit<typename std::common_type<Rep1, Rep2>::type, Scale>
    >::type
operator%(const unit<Rep1, Scale>& u, const Rep2& s)
{
    typedef typename std::common_type<Rep1, Rep2>::type Cr;
    typedef unit<Cr, Scale> Cu;
    return Cu(Cu(u).count() % static_cast<Cr>(s));
}

template <class Rep1, class Scale1, class Rep2, class Scale2>
inline constexpr typename std::common_type<unit<Rep1, Scale1>, unit<Rep2, Scale2> >::type
operator%(const unit<Rep1, Scale1>& lhs, const unit<Rep2, Scale2>& rhs)
{
    typedef typename std::common_type<Rep1, Rep2>::type Cr;
    typedef typename std::common_type<unit<Rep1, Scale1>, unit<Rep2, Scale2> >::type Cu;
    return Cu(static_cast<Cr>(Cu(lhs).count()) % static_cast<Cr>(Cu(rhs).count()));
}

} // namespace datasize

#endif //UNIT_DATASIZE_HPP
