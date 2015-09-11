set(output_files
  EOSTable.cpp
)
add_library(output ${output_files})
target_link_libraries(output eos)

