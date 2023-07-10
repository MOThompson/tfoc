# tfoc
Thin film optical calculator

This code implements the matrix method for calculating the reflectance and transmission (and absorption) of an arbitrary optical stack.
It is built primarily for the Windows environment using Visual Studio compilers, but the code should be sufficienty clean to enable porting to Linux
with minimal effort.

# Database information
The database.nk folder contains a set of n,k files that were extracted from numerous sources.  Care must be exercised to ensure that any copyright
restrictions are observed.

## Location of the database
This is a perpetual problem.  In Windows, the program looks for the database files in the Local Apps directory.
