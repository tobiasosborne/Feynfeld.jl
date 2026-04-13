set pagination off
set print elements 200
set print repeats 0

# Break inside qgen at the rule iteration moment (line 13889 - vfo lookup)
break qgraf-grind.f08:13889
commands
  silent
  printf ">>> qgen vfo lookup: vind=%d vv=%d vdeg=%d rdeg=%d\n", vind, vv, vdeg[vv-1], rdeg[vv-1]
  printf "    pmap[vv,1..vdeg]: "
  set $i = 0
  while $i < vdeg[vv-1]
    printf "%d ", pmap[vv-1+($i)*128]
    set $i = $i + 1
  end
  printf "\n"
  continue
end

# Break at the moment a rule is fully accepted (line 13900: go to 163 after rdeg loop)
break qgraf-grind.f08:13900
commands
  silent
  printf "    >> rule MATCH at vind=%d vv=%d  vfo=%d\n", vind, vv, vfo[vv-1]
  printf "       rule[1..vdeg]: "
  set $i = 1
  while $i <= vdeg[vv-1]
    printf "%d ", stib[vfo[vv-1]+$i-1]
    set $i = $i + 1
  end
  printf "\n"
  continue
end

# Break at emission point (line 13988)
break qgraf-grind.f08:13988
commands
  silent
  printf "    *** EMISSION at vind=%d (n=%d)\n", vind, n
  printf "    pmap snapshot:\n"
  set $vv2 = 1
  while $vv2 <= n
    printf "      v%d (deg %d): ", $vv2, vdeg[$vv2-1]
    set $i = 1
    while $i <= vdeg[$vv2-1]
      printf "%d ", pmap[$vv2-1+($i-1)*128]
      set $i = $i + 1
    end
    printf "\n"
    set $vv2 = $vv2 + 1
  end
  continue
end

run ctrl.dat
quit
