################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += $(wildcard ../test/*.cpp)

OBJS += $(patsubst ../test/%.cpp, ./test/%.o, $(wildcard ../test/*.cpp))

CPP_DEPS += $(patsubst ../test/%.cpp, ./test/%.d, $(wildcard ../test/*.cpp))


# Each subdirectory must supply rules for building sources it contributes
test/%.o: ../test/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/local/include -g -Wall -c -fmessage-length=0 -mtune=native -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


