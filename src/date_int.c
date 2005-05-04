#include "date.h"

#define isleap(y) \
 ((y)%4 == 0 && (y)%100 != 0 || (y)%400 == 0)

static int Dtab [2][13] =
{
  {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
  {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
};

Date *date_interval(const Date *d1, const Date *d2)
{
   static Date result;
   int months, days, years, prev_month;

   /* Compute the interval - assume d1 precedes d2 */
   years = d2->year - d1->year;
   months = d2->month - d1->month;
   days = d2->day - d1->day;

   /* Do obvious corrections (days before months!)
    *
    * This is a loop in case the previous month is
    * February, and days < -28.
    */
   prev_month = d2->month - 1;
   while (days < 0)
   {
      /* Borrow from the previous month */
      if (prev_month == 0)
         prev_month = 12;
      --months;
      days += Dtab[isleap(d2->year)][prev_month--];
   }

   if (months < 0)
   {
      /* Borrow from the previous year */
      --years;
      months += 12;
   }

   /* Prepare output */
   result.month = months;
   result.day = days;
   result.year = years;
   return &result;
}
