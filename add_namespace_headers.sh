#!/usr/bin/env bash
#
# Inserts the tadah::models namespace in header files after:
#   1) the last #include (if found and below the #define line),
#   2) otherwise after the #define line,
# and places the closing brace block exactly one line above the final #endif.
# If guard patterns are missing, fallback inserts occur at the file top and bottom.

NAMESPACE_OPEN='namespace tadah {\nnamespace models {\n'
NAMESPACE_CLOSE='} // namespace models\n} // namespace tadah\n'
FALLBACK_LOG='namespace_insertion_headers_fallback.log'

# Empty or create the fallback log
> "$FALLBACK_LOG"

if [ "$#" -lt 1 ]; then
  echo "Usage: $0 file1.h file2.h ..."
  exit 1
fi

for file in "$@"; do
  # Only process .h files
  if [[ "$file" != *.h ]]; then
    echo "Skipping $file (not .h)"
    continue
  fi

  if [ ! -f "$file" ]; then
    echo "File not found: $file"
    continue
  fi

  # Locate header guard lines
  GUARD_START_LINE=$(grep -En '^[[:space:]]*#ifndef[[:space:]]+[A-Za-z0-9_]+_H|^[[:space:]]*#ifndef[[:space:]]+[A-Za-z0-9_]+_h' "$file" \
    | head -n1 | cut -d: -f1)
  GUARD_DEFINE_LINE=$(grep -En '^[[:space:]]*#define[[:space:]]+[A-Za-z0-9_]+_H|^[[:space:]]*#define[[:space:]]+[A-Za-z0-9_]+_h' "$file" \
    | head -n1 | cut -d: -f1)

  # Locate final #endif
  ENDIF_LINE=$(grep -En '^[[:space:]]*#endif' "$file" | tail -n1 | cut -d: -f1)

  # Locate the last #include line
  LAST_INCLUDE_LINE=$(grep -En '^[[:space:]]*#include[[:space:]]*[<"]' "$file" | tail -n1 | cut -d: -f1)

  # Check for valid guard lines in ascending order
  if [ -n "$GUARD_START_LINE" ] && \
     [ -n "$GUARD_DEFINE_LINE" ] && \
     [ -n "$ENDIF_LINE" ] && \
     [ "$GUARD_START_LINE" -lt "$GUARD_DEFINE_LINE" ] && \
     [ "$GUARD_DEFINE_LINE" -lt "$ENDIF_LINE" ]; then

    # Determine insertion line for the opening namespace
    # If a last #include is found and is below #define, insert after it.
    # Otherwise, insert after the #define.
    if [ -n "$LAST_INCLUDE_LINE" ] && [ "$LAST_INCLUDE_LINE" -gt "$GUARD_DEFINE_LINE" ]; then
      INSERT_OPEN_LINE=$LAST_INCLUDE_LINE
    else
      INSERT_OPEN_LINE=$GUARD_DEFINE_LINE
    fi

    # Insert the namespace open one line after the chosen line
    INSERT_OPEN_LINE=$((INSERT_OPEN_LINE + 1))
    sed -i "${INSERT_OPEN_LINE}a\\
${NAMESPACE_OPEN}
" "$file"

    # Insert the namespace close one line above the final #endif
    sed -i "${ENDIF_LINE}i\\
${NAMESPACE_CLOSE}
" "$file"

  else
    # Fallback if guard patterns are missing or invalid
    echo "Fallback insertion in $file" | tee -a "$FALLBACK_LOG"
    # Insert at file top
    sed -i "1i\\
${NAMESPACE_OPEN}
" "$file"
    # Insert at file end
    echo -e "\n${NAMESPACE_CLOSE}" >> "$file"
  fi
done

