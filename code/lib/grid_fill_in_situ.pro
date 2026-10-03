pro grid_fill_in_situ,nyg,nz,dat

for m = 0, 11 do for z = 0, nz-1 do for y = 0, nyg-1 do for j = 0, 1 do begin
  if ~finite(dat[j,y,z,m]) then begin
    if m eq 0 and finite(dat[j,y,z,-1]) and finite(dat[j,y,z,1]) then dat[j,y,z,m] = 0.5 * (dat[j,y,z,-1] + dat[j,y,z,1])
    if m eq 11 and finite(dat[j,y,z,10]) and finite(dat[j,y,z,0]) then dat[j,y,z,m] = 0.5 * (dat[j,y,z,0] + dat[j,y,z,11])
    if m ge 1 and m le 10 then if finite(dat[j,y,z,m-1]) and finite(dat[j,y,z,m+1]) then dat[j,y,z,m] = 0.5 * (dat[j,y,z,m-1] + dat[j,y,z,m+1])
  endif
endfor
for m = 0, 11 do for z = 0, nz-1 do for y = 1, nyg-2 do for j = 0, 1 do if ~finite(dat[j,y,z,m]) then if finite(dat[j,y-1,z,m]) and finite(dat[j,y+1,z,m]) then dat[j,y,z,m] = $
  0.5 * (dat[j,y-1,z,m] + dat[j,y+1,z,m])
for m = 0, 11 do for z = 1, nz-2 do for y = 0, nyg-1 do for j = 0, 1 do if ~finite(dat[j,y,z,m]) then if finite(dat[j,y,z-1,m]) and finite(dat[j,y,z+1,m]) then dat[j,y,z,m] = $
  0.5 * (dat[j,y,z-1,m] + dat[j,y,z+1,m])
for m = 0, 11 do for z = 1, nz-2 do for y = 0, nyg-1 do for j = 0, 1 do begin
  gd = where(finite(dat[j,y,z-1:z+1,m]),ngd)
  if ngd ge 2 then dat[j,y,z,m] = mean(dat[j,y,z-1:z+1,m],/nan)
endfor
for m = 0, 11 do for z = 0, nz-1 do for y = 1, nyg-2 do for j = 0, 1 do begin
  gd = where(finite(dat[j,y-1:y+1,z,m]),ngd)
  if ngd ge 2 then dat[j,y,z,m] = mean(dat[j,y-1:y+1,z,m],/nan)
endfor

end