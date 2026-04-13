set pagination off
set print elements 200
set print repeats 0
file qgraf-grind

# Break right after qrvi finishes building dpntro and rotvpo
# qrvi sets nrot[degree] and rotvpo[degree]; we want to print these after the build
# The build completes around line 22243

break qgraf-grind.f08:22247
commands
  silent
  printf "===== After qrvi: dpntro built =====\n"
  printf "mrho=%d nrho=%d\n", mrho, nrho
  print nrot
  printf "rotvpo[3] (start of 3-vertex rule pool):\n"
  print rotvpo
  # Print the actual rule pool: at rotvpo[3] there are nrot(3)*(3+1) ints
  printf "First 20 entries of stib starting at rotvpo[3]:\n"
  set $i = 0
  while $i < 60
    printf "stib[rotvpo(3)+%d] = %d\n", $i, stib[rotvpo[3]+$i-1]
    set $i = $i + 1
  end
  continue
end

run ctrl.dat
quit
