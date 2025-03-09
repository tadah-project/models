#!/usr/bin/env bash
#
# Inserts the tadah::models namespace in .cpp source files such that:
#   - The opening namespace is immediately after the last #include line.
#   - The closing namespace is appended at the end of the file.
#
# If no #include is found, insert the namespace at line 1 and log the file to a fallback log.

NAMESPACE_OPEN="namespace tadah {\nnamespace models {\n"
NAMESPACE_CLOSE="} // namespace models\n} // namespace tadah\n"
FALLBACK_LOG="namespace_insertion_sources_fallback.log"

> "$FALLBACK_LOG"

if [ "$#" -lt 1 ]; then
  echo "Usage: $0 file1.cpp file2.cpp ..."
  exit 1
fi

for file in "$@"; do
  # Basic check for .cpp extension (feel free to allow .cc, .cxx, .C, etc.)
  if [[ "$file" != *.cpp ]]; then
    echo "Skipping $file (not a .cpp file)."
    continue
  fi

  if [ ! -f "$file" ]; then
    echo "Skipping $file (file not found)."
    continue
  fi

  # Find the last #include line
  # This regex matches lines that begin (with optional whitespace) followed by '#include' 
  # and then either <...> or "...". Adjust as you see fit.
  LAST_INCLUDE_LINE=$(grep -En '^[[:space:]]*#include[[:space:]]*[<"].+[>"]' "$file" | tail -n 1 | cut -d: -f1)

  if [ -n "$LAST_INCLUDE_LINE" ]; then
    # We found at least one #include. Insert namespace immediately after it.
    # + 1 because we want to insert *below* that line.
    INSERT_LINE=$((LAST_INCLUDE_LINE + 1))

    sed -i "${INSERT_LINE}i\\
${NAMESPACE_OPEN}
" "$file"
  else
    # Fallback: no includes found. Insert at the top of the file.
    echo "No #include found in $file. Using fallback insertion at line 1." | tee -a "$FALLBACK_LOG"
    sed -i "1i\\
${NAMESPACE_OPEN}
" "$file"
  fi

  # Insert the closing namespace at the very end
  echo -e "\n${NAMESPACE_CLOSE}" >> "$file"

done

echo "Done. Check '$FALLBACK_LOG' for any fallback insertions."

