<b>a) Która z trzech zbadanych metod okazała się najszybsza?</b><br>
Najszybsza okazała się metoda Gaussa-Seidla, dla której czas nie przekroczył dla najwyższej wartości N (12000) 0.1s. 
Niewiele gorsza okazała się metoda Jacobiego, dla której czas również nie przekroczył 0.12s.
Dla metody najgorszej - metody bezpośredniego rozwiązania równania macierzowego czas wyniósł prawie 5s, 
czyli około 40 razy więcej niż poprzednie metody.

<b>b) Która metoda iteracyjna potrzebuje mniejszej liczby iteracji, żeby zbiec się do prawidłowego wyniku?</b><br>
Mniejszej liczby iteracji w celu zbiegnięcia się do prawidłowego wyniku potrzebuje metoda metoda Gaussa-Seidla, 
dla której liczba iteracji dla najwyższej wartości N wynosiła 112. Liczba iteracji dla metody Jacobiego 
wynosiła dla tych samych danych około 2 razy więcej, a dla najwyższej wartości N 216.
<br>
<br>
<b>a) Jaka jest norma błędu rezydualnego dla każdego sposobu rozwiązania równania macierzowego?</b><br>
```norma metoda Gaussa: 6.9177e-13
norma metoda Jacobiego: 1.268575574190041e+308
minimum lokalne: 254.2894, dla liczby stron wynoszącej: 2
norma metoda Gaussa-Seidla: 64367347963.3897
minimum lokalne: 4.5044, dla liczby stron wynoszącej: 40
```

Wartość normy błędu rezydualnego przy obliczeniach metodą Jacobiego <br>
po pierwszej iteracji spadła osiągając minimum lokalne.
Natomiast przy obliczeniach metodą Gaussa-Seidla początkowo znacząco zmniejszała się
lecz od ok 130 iteracji zaczęła wzrastać.<br>

<b>b) Czy metody iteracyjne zbiegają się?</b><br>
Dla problemu analizy elektromagnetycznej filtru mikrofalowego metody iteracyjne <br>
Jacobiego oraz Gaussa-Seidla nie zbiegają się.<br>
