#include "BEMUnit.h"
#include "InputUnit.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

int main()
{
  namespace fs = std::filesystem;

  const fs::path sourceInput = fs::path(MKBEM_TEST_DATA_DIR) / "THM_TeST2.dat";
  const fs::path workInput = fs::current_path() / "THM_TeST2.dat";
  const fs::path workOutput = fs::current_path() / "THM_TeST2.out";

  std::error_code error;
  fs::remove(workOutput, error);
  fs::copy_file(sourceInput, workInput, fs::copy_options::overwrite_existing, error);
  if (error)
  {
    std::cerr << "Failed to prepare benchmark input: " << error.message() << std::endl;
    return 1;
  }

  TBEMInput input;
  if (!input.ReadFile(workInput.string()))
  {
    std::cerr << "Failed to parse benchmark input" << std::endl;
    return 1;
  }

  const int expectedTimeSteps = input.GetData().NT;
  try
  {
    TTHMBEM model(input);
    model.Execute();
  }
  catch (const std::exception &error)
  {
    std::cerr << "Benchmark execution failed: " << error.what() << std::endl;
    return 1;
  }

  std::ifstream output(workOutput);
  if (!output.is_open())
  {
    std::cerr << "Benchmark output was not created: " << workOutput << std::endl;
    return 1;
  }

  std::ostringstream buffer;
  buffer << output.rdbuf();
  std::string text = buffer.str();
  std::transform(text.begin(), text.end(), text.begin(), [](unsigned char ch) {
    return static_cast<char>(std::tolower(ch));
  });

  if (text.find("nan") != std::string::npos || text.find("inf") != std::string::npos)
  {
    std::cerr << "Benchmark output contains a non-finite value" << std::endl;
    return 1;
  }

  int timeStepCount = 0;
  std::string::size_type pos = 0;
  while ((pos = text.find("time =", pos)) != std::string::npos)
  {
    timeStepCount++;
    pos += 6;
  }

  if (timeStepCount != expectedTimeSteps)
  {
    std::cerr << "Expected " << expectedTimeSteps << " output time steps, found " << timeStepCount << std::endl;
    return 1;
  }

  return 0;
}