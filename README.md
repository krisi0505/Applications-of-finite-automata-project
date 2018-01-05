# Applications of finite automata project
The goal of the project was to construct a minimal deterministic finite automaton by a given regular expression extended with operations for intersection, complement, reversal and difference. 

I have used the power set construction algorithm for determinisation and the Hopcroft's algorithm for minimization. I managed to do the extended operations using De Morgan's laws.

The example automata accepts only valid and existing dates i.e. it doesn't accept dates such as 'FEBRUARY 29, 2017$'. The dates are accepted in the format 'MONTH DATE, YEAR$', where:
* MONTH is JANUARY, FEBRUARY, MARCH, APRIL, MAY, JUNE, JULY, AUGUST, SEPTEMBER, OCTOBER, NOVEMBER or DECEMBER;
* DATE is a number between 1 and 31;
* YEAR is a year between 1 and positive infinity
