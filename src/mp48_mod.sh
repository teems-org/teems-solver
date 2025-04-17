#!/bin/bash
cp ./hsl_mp48d.f90 ./hsl_mp48d_original.f90
cp ./ddeps.f ./ddeps_original.f
sed -i 's/! COPYRIGHT (c) 2003 Council for the Central Laboratory/!-Add MYIW1 and MYLMX;'"\n"'!-Change CGLOB=Max col number in a block (including border) x NBLOCK'"\n"'! COPYRIGHT (c) 2003 Council for the Central Laboratory/g' ./hsl_mp48d.f90
sed -i 's/      INTEGER RINTER         ! No. of rows in interface problem/      INTEGER RINTER         ! No. of rows in interface problem'"\n"'      INTEGER MAXSBCOLS      ! Maximum No. of cols in single block (including border part)/g' ./hsl_mp48d.f90
sed -i 's/         ALLOCATE (data%private%CGLOB(NEQ,NBLOCK), \&/         ALLOCATE (data%private%CGLOB(data%MAXSBCOLS,NBLOCK), \& !NEQ/g' ./hsl_mp48d.f90
sed -i 's/      INTEGER IOS        !  Holds IOSTAT parameter for files./      INTEGER(8) MYLMX'"\n"'      INTEGER IOS        !  Holds IOSTAT parameter for files./g' ./hsl_mp48d.f90

sed -i '0,/      INTEGER, DIMENSION (:), ALLOCATABLE :: IW1/s/      INTEGER, DIMENSION (:), ALLOCATABLE :: IW1/      INTEGER, DIMENSION (:), ALLOCATABLE :: IW1'"\n"'      INTEGER(8), DIMENSION (:), ALLOCATABLE :: MYIW1/g' ./hsl_mp48d.f90

sed -i '/                ROWPTR(1:LORDER+1), \&/{ N; s/                ROWPTR(1:LORDER+1), \&\n                STAT=ST)/                ROWPTR(1:LORDER+1), \&\n                STAT=ST)\n      DEALLOCATE (MYIW1,STAT=ST)\n      ALLOCATE (MYIW1(1:NEQ),STAT=ST)/g }' ./hsl_mp48d.f90
sed -i '/      IW1(1:NEQ) = 0/{ N; s/      IW1(1:NEQ) = 0\n      LMX = 0/      MYIW1(1:NEQ) = 0\n      MYLMX = 0/g }' ./hsl_mp48d.f90
sed -i '/         IDUP1 = 0/{ N; s/         IDUP1 = 0\n         IW1(1:NEQ) = 0/         IDUP1 = 0\n         MYIW1(1:NEQ) = 0/g }' ./hsl_mp48d.f90
sed -i '/             J = J + 1/{ N; s/             J = J + 1\n             IF (IW1(I) < L+LMX) THEN/             J = J + 1\n             IF (MYIW1(I) < L+MYLMX) THEN/g }' ./hsl_mp48d.f90
sed -i '/               JJ = JJ + 1/{ N; s/               JJ = JJ + 1\n               IW1(I) = L + LMX/               JJ = JJ + 1\n               MYIW1(I) = L + MYLMX/g }' ./hsl_mp48d.f90
sed -i '/  230    CONTINUE/{ N; s/  230    CONTINUE\n         LMX = LMX + NEQ/  230    CONTINUE\n         MYLMX = MYLMX + NEQ/g }' ./hsl_mp48d.f90
sed -i '/      IF (data%IDUP > 0) IW1(1:NEQ) = 0/{ N; s/      IF (data%IDUP > 0) IW1(1:NEQ) = 0\n      LMX = 0/      IF (data%IDUP > 0) MYIW1(1:NEQ) = 0\n      MYLMX = 0/g }' ./hsl_mp48d.f90

sed -i '0,/                 IF (IW1(I) < L+LMX) THEN/{ N; s/                 IF (IW1(I) < L+LMX) THEN\n                    IW1(I) = L + LMX/                 IF (MYIW1(I) < L+MYLMX) THEN\n                    MYIW1(I) = L + MYLMX/g }' ./hsl_mp48d.f90

sed -i '/  270        CONTINUE/{ N; s/  270        CONTINUE\n             LMX = LMX + NEQ/  270        CONTINUE\n             MYLMX = MYLMX + NEQ/g }' ./hsl_mp48d.f90
sed -i '/  271        CONTINUE/{ N; s/  271        CONTINUE\n             LMX = LMX + NEQ/  271        CONTINUE\n             MYLMX = MYLMX + NEQ/g }' ./hsl_mp48d.f90
sed -i '/      DEALLOCATE (IW1,STAT=ST)/{ N; s/      DEALLOCATE (IW1,STAT=ST)\n      DEALLOCATE (IRN,STAT=ST)/      DEALLOCATE (IW1,STAT=ST)\n      DEALLOCATE (MYIW1,STAT=ST)\n      DEALLOCATE (IRN,STAT=ST)/g }' ./hsl_mp48d.f90

sed -i '/                 IF (IW1(I) < L+LMX) THEN/{ N; s/                 IF (IW1(I) < L+LMX) THEN\n                    IW1(I) = L + LMX/                 IF (MYIW1(I) < L+MYLMX) THEN\n                    MYIW1(I) = L + MYLMX/g }' ./hsl_mp48d.f90
sed -i '0,/                 IF (MYIW1(I) < L+MYLMX) THEN/{ N; s/                 IF (MYIW1(I) < L+MYLMX) THEN\n                    MYIW1(I) = L + MYLMX/                 IF (IW1(I) < L+LMX) THEN\n                    IW1(I) = L + LMX/g }' ./hsl_mp48d.f90
sed -i '0,/                 IF (MYIW1(I) < L+MYLMX) THEN/{ N; s/                 IF (MYIW1(I) < L+MYLMX) THEN\n                    MYIW1(I) = L + MYLMX/                 IF (IW1(I) < L+LMX) THEN\n                    IW1(I) = L + LMX/g }' ./hsl_mp48d.f90
sed -i '3394,3437 s/MYIW1/IW1/g; 3394,3437 s/MYLMX/LMX/g' ./hsl_mp48d.f90
sed -i 's/MA48DD/ZA48DD/g' ./ddeps.f
sed -i 's/MA48DD/ZA48DD/g' ./hsl_mp48d.f90
sed -i 's/MA50CD/ZA50CD/g' ./ddeps.f
sed -i 's/MA50CD/ZA50CD/g' ./hsl_mp48d.f90
sed -i 's/MA50HD/ZA50HD/g' ./ddeps.f
sed -i 's/MA50HD/ZA50HD/g' ./hsl_mp48d.f90
sed -i 's/MC71AD/ZC71AD/g' ./ddeps.f
sed -i 's/MC71AD/ZC71AD/g' ./hsl_mp48d.f90
sed -i 's/MA50BD/ZA50BD/g' ./ddeps.f
sed -i 's/MA50BD/ZA50BD/g' ./hsl_mp48d.f90
sed -i 's/MA50GD/ZA50GD/g' ./ddeps.f
sed -i 's/MA50GD/ZA50GD/g' ./hsl_mp48d.f90
sed -i 's/MA50FD/ZA50FD/g' ./ddeps.f
sed -i 's/MA50FD/ZA50FD/g' ./hsl_mp48d.f90
sed -i 's/MA50ED/ZA50ED/g' ./ddeps.f
sed -i 's/MA50ED/ZA50ED/g' ./hsl_mp48d.f90
sed -i 's/MA50AD/ZA50AD/g' ./ddeps.f
sed -i 's/MA50AD/ZA50AD/g' ./hsl_mp48d.f90
sed -i 's/MA50DD/ZA50DD/g' ./ddeps.f
sed -i 's/MA50DD/ZA50DD/g' ./hsl_mp48d.f90
sed -i 's/MC13DD/ZC13DD/g' ./ddeps.f
sed -i 's/MC13DD/ZC13DD/g' ./hsl_mp48d.f90
sed -i 's/MC13ED/ZC13ED/g' ./ddeps.f
sed -i 's/MC13ED/ZC13ED/g' ./hsl_mp48d.f90
sed -i 's/MC21AD/ZC21AD/g' ./ddeps.f
sed -i 's/MC21AD/ZC21AD/g' ./hsl_mp48d.f90
sed -i 's/MC21BD/ZC21BD/g' ./ddeps.f
sed -i 's/MC21BD/ZC21BD/g' ./hsl_mp48d.f90
sed -i 's/MA48AD/ZA48AD/g' ./ddeps.f
sed -i 's/MA48AD/ZA48AD/g' ./hsl_mp48d.f90
sed -i 's/MA48BD/ZA48BD/g' ./ddeps.f
sed -i 's/MA48BD/ZA48BD/g' ./hsl_mp48d.f90
sed -i 's/MA48CD/ZA48CD/g' ./ddeps.f
sed -i 's/MA48CD/ZA48CD/g' ./hsl_mp48d.f90
sed -i 's/MA48ID/ZA48ID/g' ./ddeps.f
sed -i 's/MA48ID/ZA48ID/g' ./hsl_mp48d.f90
sed -i 's/MA50ID/ZA50ID/g' ./ddeps.f
sed -i 's/MA50ID/ZA50ID/g' ./hsl_mp48d.f90
sed -i 's/MC29AD/ZC29AD/g' ./ddeps.f
sed -i 's/MC29AD/ZC29AD/g' ./hsl_mp48d.f90
sed -i 's/MC59DD/ZC59DD/g' ./ddeps.f
sed -i 's/MC59DD/ZC59DD/g' ./hsl_mp48d.f90
sed -i 's/MC59FD/ZC59FD/g' ./ddeps.f
sed -i 's/MC59FD/ZC59FD/g' ./hsl_mp48d.f90
sed -i 's/MC59CD/ZC59CD/g' ./ddeps.f
sed -i 's/MC59CD/ZC59CD/g' ./hsl_mp48d.f90
sed -i 's/MC59ED/ZC59ED/g' ./ddeps.f
sed -i 's/MC59ED/ZC59ED/g' ./hsl_mp48d.f90
sed -i 's/MC59BD/ZC59BD/g' ./ddeps.f
sed -i 's/MC59BD/ZC59BD/g' ./hsl_mp48d.f90
sed -i 's/MC59AD/ZC59AD/g' ./ddeps.f
sed -i 's/MC59AD/ZC59AD/g' ./hsl_mp48d.f90