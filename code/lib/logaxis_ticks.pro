pro logaxis_ticks,p1,p2,tickvals

;  Procedure to make tick labels for a plot with a log axis, typically pressure.

n1 = fix((p1+1)/100)
if n1 lt 0 then n1 = 0
n2 = fix(10 - p2/10)
if n2 lt 0 then n2 = 0
if n2 eq 10 then n2-=1
n3 = fix(10 - p2)
if n3 le 0 then n3 = 0
nticks = n1 + n2 + n3
tickvals = intarr(nticks)

if n1 gt 0 then tickvals[0:n1-1] = indgen(n1)*(-100) + n1*100
if n2 gt 0 then tickvals[n1:n1+n2-1] = indgen(n2)*(-10) + 90
if n3 gt 0 then tickvals[n1+n2:nticks-1] = indgen(n3)*(-1) + 9

test = where(tickvals lt p2,ntest)
if ntest gt 0 then tickvals = tickvals[where(tickvals ge p2)]

; ----  Add a tick for 150 hPa.
if p1 ge 200 and p2 le 100 then begin
	nticks = nticks + 1
	tickvals2 = intarr(nticks)
	i1 = where(tickvals ge 200,ni1)
	tickvals2[i1] = tickvals[i1]
	tickvals2[ni1] = 150
	i2 = where(tickvals le 100,ni2)
	tickvals2[i2+1] = tickvals[i2]
	tickvals = tickvals2
endif

end
