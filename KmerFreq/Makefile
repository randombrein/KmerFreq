################################################################################
## Program that produce the most frequent DNA-k-mers of arbitrary length sorted
## by frequency in a given FASTQ file.
################################################################################
TARGET=kmer_freq
SRC_DIR=src
INC_DIR=include
BIN_DIR=bin

SRCS=$(wildcard $(SRC_DIR)/*.cpp)
OBJS=$(patsubst %, $(BIN_DIR)/%, $(notdir $(SRCS:.cpp=.o)))

CC=gcc
CXX=g++
RM=rm -f
MKDIR_P = mkdir -p
CPPFLAGS=-Wall -std=gnu++11 -O3 -I$(INC_DIR)
LDFLAGS=$(shell root-config --ldflags)
LDLIBS=$(shell root-config --libs)

GREEN=\033[1;32m
YELLOW=\033[1;33m
NC=\033[1;0m


$(BIN_DIR)/$(TARGET) : $(OBJS)
	@echo "$(YELLOW)linking...$(NC)"
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)
	@echo "$(GREEN)finished!$(NC)"

$(BIN_DIR)/%.o : $(SRC_DIR)/%.cpp
	@echo "$(YELLOW)compiling...$(NC)"
	$(CXX) $(CPPFLAGS) -c $< -o $@

depend: .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;

.PHONY : help clean run

run : $(BIN_DIR)/$(TARGET)
	@echo "$(GREEN)running...$(NC)"
	$(BIN_DIR)/$(TARGET) --filename res/ERR037900_1.first1000.fastq --kmersize 30 --topcount 25


${BIN_DIR}:
	$(MKDIR_P) $(BIN_DIR)

clean :
	$(RM) $(OBJS) $(BIN)/$(TARGET)

dist-clean : clean
	$(RM) *~ .depend

include .depend
