################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/lifFile.cpp \
../src/locatorfromlif.cpp \
../src/multiscalefinder.cpp \
../src/octavefinder.cpp \
../src/reconstructor.cpp \
../src/traj.cpp 

OBJS += \
./src/lifFile.o \
./src/locatorfromlif.o \
./src/multiscalefinder.o \
./src/octavefinder.o \
./src/reconstructor.o \
./src/traj.o 

CPP_DEPS += \
./src/lifFile.d \
./src/locatorfromlif.d \
./src/multiscalefinder.d \
./src/octavefinder.d \
./src/reconstructor.d \
./src/traj.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/local/include -O3 -g3 -Wall -c -fopenmp -fmessage-length=0 -mtune=native -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


