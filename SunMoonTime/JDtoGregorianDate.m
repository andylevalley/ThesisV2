function UTC = JDtoGregorianDate(JD)

% Taken from Vallado, Algorithm 22

T_1900 = (JD-2415019.5)/365.25;
yr = 1900+floor(T_1900);
LeapYrs = floor((yr-1900-1)*.25);
days = (JD-2415019.5)-((yr-1900)*365+LeapYrs);
if days<1
    yr = yr-1;
    LeapYrs = floor((yr-1900-1)*.25);
    days = (JD-2415019.5)-((yr-1900)*365+LeapYrs);
end
LMonth = [31 28 31 30 31 30 31 31 30 31 30 31];
if mod(yr,4)==0
    LMonth(2) = 29;
end
dayofyr = floor(days);
mo_counter = 1;
summer = 0;
while (dayofyr > (summer + LMonth(mo_counter))) && (mo_counter < 12)
    summer = summer + LMonth(mo_counter);
    mo_counter = mo_counter + 1;
end
mo = mo_counter;
day = dayofyr-summer;
tau = (days-dayofyr)*24;
h = floor(tau);
min = floor((tau-h)*60);
s = (tau-h-min/60)*3600;

UTC(1) = yr;
UTC(2) = mo;
UTC(3) = day;
UTC(4) = h;
UTC(5) = min;
UTC(6) = s;