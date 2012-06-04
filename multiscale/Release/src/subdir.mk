################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += $(wildcard ../src/*.cpp) 

OBJS += $(patsubst ../src/%.cpp, ./src/%.o, $(wildcard ../src/*.cpp))

CPP_DEPS += $(patsubst ../src/%.cpp, ./src/%.d, $(wildcard ../src/*.cpp))


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/local/include -O3 -g3 -Wall -c -fopenmp -fmessage-length=0 -mtune=native -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


