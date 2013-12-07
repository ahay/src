	IMPLICIT NONE
	INTEGER, DIMENSION(1) :: OLD, SEED ! THIS PROGRAM ASSUMES K = 1
	INTEGER :: I, K
	REAL, DIMENSION(3) :: HARVEST
	SEED(1) = 12345
	CALL RANDOM_SEED
	CALL RANDOM_SEED(SIZE=K)
	WRITE(*,*) ' Number of integers for starting value = ', K
	CALL RANDOM_SEED(GET=OLD(1:K))
	WRITE(*,*) ' Old starting value = ', OLD
	CALL RANDOM_NUMBER(HARVEST)
	WRITE(*,*) ' Random numbers : ', HARVEST
	CALL RANDOM_SEED(GET=OLD(1:K))
	WRITE(*,*) ' Present starting value = ', OLD
	CALL RANDOM_SEED(PUT=SEED(1:K))
	CALL RANDOM_SEED(GET=OLD(1:K))
	WRITE(*,*) ' New starting value = ', OLD
	CALL RANDOM_NUMBER(HARVEST)
	WRITE(*,*) ' Random numbers : ', HARVEST
	DO I = 1, 3
		CALL RANDOM_SEED(GET=OLD(1:K))
		WRITE(*,*) ' Present starting value = ', OLD
		CALL RANDOM_NUMBER(HARVEST)
		WRITE(*,*) ' Random numbers : ', HARVEST
		CALL RANDOM_NUMBER(HARVEST)
		WRITE(*,*) ' Random numbers : ', HARVEST
	END DO
	END
