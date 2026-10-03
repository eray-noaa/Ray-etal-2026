pro grid_fill_in_situ_seas,nyg,nz,dat

for m = 0, 3 do for z = 0, nz-1 do for y = 1, nyg-3 do for j = 0, 1 do if ~finite(dat[j,y,z,m]) and ~finite(dat[j,y+1,z,m]) then $
  if finite(dat[j,y-1,z,m]) and finite(dat[j,y+2,z,m]) then dat[j,y,z,m] = 0.67 * dat[j,y-1,z,m] + 0.33 * dat[j,y+2,z,m]
for m = 0, 3 do for z = 0, nz-1 do for y = 1, nyg-2 do for j = 0, 1 do if ~finite(dat[j,y,z,m]) then if finite(dat[j,y-1,z,m]) and finite(dat[j,y+1,z,m]) then dat[j,y,z,m] = $
  0.5 * (dat[j,y-1,z,m] + dat[j,y+1,z,m])
for m = 0, 3 do for z = 1, nz-3 do for y = 0, nyg-1 do for j = 0, 1 do if ~finite(dat[j,y,z,m]) and ~finite(dat[j,y,z+1,m]) then $
  if finite(dat[j,y,z-1,m]) and finite(dat[j,y,z+2,m]) then dat[j,y,z,m] = 0.67 * dat[j,y,z-1,m] + 0.33 * dat[j,y,z+2,m]
for m = 0, 3 do for z = 1, nz-2 do for y = 0, nyg-1 do for j = 0, 1 do if ~finite(dat[j,y,z,m]) then if finite(dat[j,y,z-1,m]) and finite(dat[j,y,z+1,m]) then dat[j,y,z,m] = $
  0.5 * (dat[j,y,z-1,m] + dat[j,y,z+1,m])
for m = 0, 3 do for z = 1, nz-2 do for y = 0, nyg-1 do for j = 0, 1 do begin
  gd = where(finite(dat[j,y,z-1:z+1,m]),ngd)
  if ngd ge 2 then dat[j,y,z,m] = mean(dat[j,y,z-1:z+1,m],/nan)
endfor
for m = 0, 3 do for z = 0, nz-1 do for y = 1, nyg-2 do for j = 0, 1 do begin
  gd = where(finite(dat[j,y-1:y+1,z,m]),ngd)
  if ngd ge 2 then dat[j,y,z,m] = mean(dat[j,y-1:y+1,z,m],/nan)
endfor

for m = 0, 3 do for z = 0, nz-1 do for y = 0, nyg-1 do for j = 0, 1 do begin
  if ~finite(dat[j,y,z,m]) then begin
    if m eq 0 and finite(dat[j,y,z,-1]) and finite(dat[j,y,z,1]) then dat[j,y,z,m] = 0.5 * (dat[j,y,z,-1] + dat[j,y,z,1])
    if m eq 3 and finite(dat[j,y,z,2]) and finite(dat[j,y,z,0]) then dat[j,y,z,m] = 0.5 * (dat[j,y,z,0] + dat[j,y,z,2])
    if m eq 1 or m eq 2 then if finite(dat[j,y,z,m-1]) and finite(dat[j,y,z,m+1]) then dat[j,y,z,m] = 0.5 * (dat[j,y,z,m-1] + dat[j,y,z,m+1])
  endif
endfor

for m = 0, 3 do for z = 0, nz-1 do for y = 1, nyg-3 do for j = 0, 1 do if ~finite(dat[j,y,z,m]) and ~finite(dat[j,y+1,z,m]) then $
  if finite(dat[j,y-1,z,m]) and finite(dat[j,y+2,z,m]) then dat[j,y,z,m] = 0.67 * dat[j,y-1,z,m] + 0.33 * dat[j,y+2,z,m]
for m = 0, 3 do for z = 0, nz-1 do for y = 1, nyg-2 do for j = 0, 1 do if ~finite(dat[j,y,z,m]) then if finite(dat[j,y-1,z,m]) and finite(dat[j,y+1,z,m]) then $
  dat[j,y,z,m] = 0.5 * (dat[j,y-1,z,m] + dat[j,y+1,z,m])
for m = 0, 3 do for z = 1, nz-3 do for y = 0, nyg-1 do for j = 0, 1 do if ~finite(dat[j,y,z,m]) and ~finite(dat[j,y,z+1,m]) then $
  if finite(dat[j,y,z-1,m]) and finite(dat[j,y,z+2,m]) then dat[j,y,z,m] = 0.67 * dat[j,y,z-1,m] + 0.33 * dat[j,y,z+2,m]
for m = 0, 3 do for z = 1, nz-2 do for y = 0, nyg-1 do for j = 0, 1 do if ~finite(dat[j,y,z,m]) then if finite(dat[j,y,z-1,m]) and finite(dat[j,y,z+1,m]) then $
  dat[j,y,z,m] = 0.5 * (dat[j,y,z-1,m] + dat[j,y,z+1,m])
for m = 0, 3 do for z = 1, nz-2 do for y = 0, nyg-1 do for j = 0, 1 do begin
  gd = where(finite(dat[j,y,z-1:z+1,m]),ngd)
  if ngd ge 2 then dat[j,y,z,m] = mean(dat[j,y,z-1:z+1,m],/nan)
endfor
for m = 0, 3 do for z = 0, nz-1 do for y = 1, nyg-2 do for j = 0, 1 do begin
  gd = where(finite(dat[j,y-1:y+1,z,m]),ngd)
  if ngd ge 2 then dat[j,y,z,m] = mean(dat[j,y-1:y+1,z,m],/nan)
endfor

for m = 0, 3 do for z = 0, nz-1 do for y = 0, nyg-1 do for j = 0, 1 do begin
  if ~finite(dat[j,y,z,m]) then begin
    if m eq 0 and finite(dat[j,y,z,-1]) and finite(dat[j,y,z,1]) then dat[j,y,z,m] = 0.5 * (dat[j,y,z,-1] + dat[j,y,z,1])
    if m eq 3 and finite(dat[j,y,z,2]) and finite(dat[j,y,z,0]) then dat[j,y,z,m] = 0.5 * (dat[j,y,z,0] + dat[j,y,z,2])
    if m eq 1 or m eq 2 then if finite(dat[j,y,z,m-1]) and finite(dat[j,y,z,m+1]) then dat[j,y,z,m] = 0.5 * (dat[j,y,z,m-1] + dat[j,y,z,m+1])
  endif
endfor

for a = 0, 5 do for m = 0, 3 do for y = 0, nyg-1 do for j = 0, 1 do begin
  gd = where(finite(dat[j,y,*,m]),ngd)
  if ngd ge 6 then begin
    dz = dat[j,y,gd,m] - shift(dat[j,y,gd,m],1)
    if max(gd) le 42 then begin
      dat[j,y,gd[-1]+1,m] = dat[j,y,gd[-1],m] + mean(dz[-2:-1])
    endif
  endif
endfor

for m = 0, 3 do for z = 1, nz-2 do for y = 0, nyg-1 do for j = 0, 1 do begin
  gd = where(finite(dat[j,y,z-1:z+1,m]),ngd)
  if ngd ge 2 then dat[j,y,z,m] = mean(dat[j,y,z-1:z+1,m],/nan)
endfor
for m = 0, 3 do for z = 0, nz-1 do for y = 1, nyg-2 do for j = 0, 1 do begin
  gd = where(finite(dat[j,y-1:y+1,z,m]),ngd)
  if ngd ge 2 then dat[j,y,z,m] = mean(dat[j,y-1:y+1,z,m],/nan)
endfor

end