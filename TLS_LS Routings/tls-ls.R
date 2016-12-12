TLS_Gen <- function (mode, scaling, ef, n, m, In, y, p, X, E, S, Corr, sigma, s0, s0_dis, err)
{
  #---------------------------------------------------------------
  #     (TLS, МНК) в общем виде
  
  #integer, intent(in) :: mode, scaling, ef, n, m
  #doubleprecision, intent(in) :: In(n,m),y(n),p(n)
  #doubleprecision, intent(out):: x(m), E(n), S(m+1), Corr(m,m), sigma(m), s0, s0_dis, err;  
  
  #doubleprecision, allocatable :: a(:,:), aa(:,:), b(:,:), y1(:), U(:,:), V(:,:), dn(:,:), cn(:,:), cn_inv(:,:)
  #doubleprecision, allocatable :: b11(:,:), b12(:,:), b21(:,:), b22(:,:), b_shur(:,:), b_inv(:,:)  
  #doubleprecision, allocatable :: b12v(:), b21v(:)
  #integer i,j, nef
  #doubleprecision ee(m,m),c(m,m),d(m,m),d1(m,m),det,E1, b11s, Cov(m,m);
  #doubleprecision ss(m+1, m+1), x2(m+1), s_m(m)
  #complex(kind(1d0)) s_c(m+1) 
  #complex(kind(1d0)), allocatable :: s_b(:)
  
  if ((m>n) | (m==n) | (n<3) | (m<2))
  { 
    err<<- (-65535)
    return;
  }
  
  
  #s<-0; x<-0; e<-0; corr<-0; sigma=0; s0=0; s0_dis=0;
  
  
  
  #========= вычисление сингулярных чисел матрицы исходных данных для оценки числа обусловленности  
  allocate(a(n,m), U(n, n), V(m, m));
  a = In;
  call lin_svd(a, s_m, U, V);
  deallocate(a, U, V);
  
  
  !========= Подготовка исходных данных, вычисляем дополненную матрицу условных уравнений
  allocate(a(n,m+1), aa(n,m+1), U(n,n), V(m+1, m+1));
  
  a = reshape(In, (/n, m+1/), y);
  
  !========= Вычисляет коэффициенты нормализации
  if (scaling>0) then
  
  allocate (dn(n,n), cn(m+1, m+1), cn_inv(m+1, m+1));   
  
  ! по строчкам   
  dn = 0;
  if ((scaling==1).or.(scaling==3)) then    
  do i=1,n
  dn(i,i) = 1/nrm2(a(i,:));
  !     dn(i,i) = 1;
  enddo
  a = matmul(dn, a);
  endif 
  
  ! по столбцам
  cn = 0;
  if ((scaling==2).or.(scaling==3)) then
  do i=1,m+1
  cn(i,i) = 1/nrm2(a(:,i));
  !     cn(i,i) = 1;
  enddo
  a = matmul(a, cn);
  cn_inv = 0;
  call Inverse(cn, cn_inv, m+1);
  endif
  
  ! норма Фробениуса
  if (scaling==4) then
  !  call NR2RR (a, fn)
  endif
  
  endif
  
  
  !=====   вычислнение решения
  select case (mode)
  case (1)
  !============ Решение через сингуляргое разложение
  call lin_svd(a, s, U, V); 
  
  !----  проверка решения
  SS = 0;
  do i=1,m    ! тут нет ошибки, последняя ячейка должна быть нулевой
  SS(i,i) = s(i);
  enddo
  aa = matmul(U,matmul(SS,V));
  ee = aa-a;
  !----  вычисляем х
  x2 =  v(:,m+1);  
  
  if ((scaling==2).or.(scaling==3)) then
  x2 = matmul(cn, x2); 
  endif;
  
  x2 = x2/x2(m+1)
  x = -1*x2;
  
  !---- Вычисление числа обусловленности
  err = (s_m(1)-s(m+1))/(s_m(m)-s(m+1));
  
  !----  вычисление sigma_0, The Total Least Squares Problem, 8.15
  s0 = s(m+1) / sqrt(n+0.0);
  !----
    
    !---- вычисление ковариационной матрицы, The Total Least Squares Problem, page 242, 8.47
  d = 0;
  do i=1,m    
  d(i,i) = s(m+1)**2;
  enddo
  
  d1 = matmul(transpose(In),In) - d;
  call inverse(d1, Cov, m)
  Cov = (1+sum(x**2)) * s0**2 * Cov;
  
  !===========
    case (2)
  !===========  Решение через собственные числа
  ss=matmul(transpose(a),a);
  
  if (ef<1) then  ! вариант польного TLS, все переменные содержат ошибки
  nef=m+1;
  
  call lin_eig_gen(ss, s_c);  
  ee = 0; s=0;
  do i=1,m
  ee(i,i) = s_c(m+1);
  s(i) = s_c(i);
  enddo
  s(m+1) = s_c(m+1);
  if (s(m+1)<0) then 
  s(m+1) = 0;
  endif
  
  
  !----  вычисление sigma_0, The Total Least Squares Problem, 8.15
  s0 = sqrt(s(m+1))/sqrt(n+0.0);
  !---- Вычисление числа обусловленности
  err =  (s_m(1)-sqrt(s_c(m+1))) / (s_m(m)-sqrt(s_c(m+1)));    
  
  !---- вычисление ковариационной матрицы, по аналогии с The Total Least Squares Problem, page 242, 8.47
  d = 0;
  do i=1,m    
  d(i,i) = s(m+1);
  enddo
  
  d1 = matmul(transpose(In),In) - d;
  call inverse(d1, Cov, m)
  
  !    этот частный случай (когда только первый столбец error-free) вошел в более общий, реализованный ниже. Код рабочий, на всякий случай оставляю
  !    elseif (ef==-1) then     
  !     nef=m;
  !     allocate (b21v(m),b12v(m),b22(m, m), b_shur(m, m), s_b(m));
  !     b21v = ss(1,2:m+1);
  !     b12v = ss(2:m+1, 1);
  !     b22 = ss(2:m+1,2:m+1);
  !     do i=1,m
  !       do j=1,m
  !         b_shur(i,j) = b21v(i)*b12v(j)
  !       enddo
  !     enddo
  
  !     b_shur = b_shur/ss(1,1);
  !     b_shur = b22 - b_shur
  
  !     call lin_eig_gen(b_shur, s_b);     
  !     ee = 0;
  !     s=0;
  !     do i=2,m
  !       ee(i,i) = s_b(m);      
  !       s(i-1) = s_b(i-1)
  !     enddo   
  !     s(m) = s_b(m);
  !     deallocate(b21v, b12v, b22, b_shur, s_b);
  
  else   ! вариант TLS-LS, когда первый ef-столбцов не содержат ошибки, а последние nef - содержат
  nef = m+1-ef; 
  allocate (b21(nef,ef), b12(ef,nef),b22(nef, nef), b_shur(nef, nef),b_inv(ef,ef), s_b(nef), b11(ef,ef));
  b12 = ss(1:ef,ef+1:m+1);
  b21 = ss(ef+1:m+1, 1:ef);
  b22 = ss(ef+1:m+1,ef+1:m+1);
  b11 = ss(1:ef,1:ef);
  
  call inverse(b11, b_inv, ef)
  ! B22-B21*Inverse(B11)*B12    
  b_shur = b22 - matmul(matmul(b21,b_inv),b12);
  
  s=0;
  call lin_eig_gen(b_shur, s_b);     
  ee = 0;
  do i=ef+1,m
  ee(i,i) = s_b(nef);             
  enddo   
  
  do i=1,nef
  s(i) = s_b(i);
  enddo
  
  !----  вычисление sigma_0, по аналогии с The Total Least Squares Problem, 8.15, возможно не правильно, надо как-то проверить
  s0 = sqrt(s_b(nef))/sqrt(n+0.0);
  
  !---- вычисление ковариационной матрицы, по аналогии с The Total Least Squares Problem, page 242, 8.47, возможно не правильно, надо как-то проверить
  d = 0;
  do i=ef+1,m    
  d(i,i) = s_b(nef)**2;
  enddo
  d1 = matmul(transpose(In),In) - d;
  call inverse(d1, Cov, m)
  
  !---- Вычисление числа обусловленности
  err = 0;
  !     err =  (s_m(1)-sqrt(s(nef))) / (s_m(m)-sqrt(s(nef)));   формула не правильная, ошибка в знаменателе. Не понятно как оценивать обусловленность в таком случае..
  
  deallocate(b21, b12, b22, b11, b_shur, b_inv, s_b);      
  endif
  
  
  !--------- Вычисление матрицы нормальных уравнений
  c=matmul(transpose(a(:,1:m)),a(:,1:m))-ee;  ! случай с масштабированием
  !c=matmul(transpose(In),In)-ee;  ! это правильно, если нет мастабирования 
  
  !----------------------------------  
    !  вычисление обратной матрицы
  ee=0;
  do i=1,m
  ee(i,i)=1;
  enddo  
  call lin_sol_gen(c, ee, d) 
  ee = matmul(c,d);  ! проверка
  
  !----   вычисление решения  
  x = matmul(d, matmul(transpose(a(:,1:m)),a(:,m+1))); ! с учетом мастабирования
  !x = matmul(d, matmul(transpose(a),y));    ! а - дополненная матрица, может быть тут все таки должен быть In? Хотя работает правильно
  
  if ((scaling==2).or.(scaling==3)) then
  x = matmul(cn, x)/cn(m+1,m+1); 
  endif;
  
  !---- заканчиваем вычисление корреляции
  Cov = (1+sum(x**2)) * s0**2 * Cov;
  
  !========
    endselect
  
  deallocate(a, aa, u, v);
  
  if (scaling>0) then
  deallocate(dn, cn, cn_inv);
  endif
  
  
  !----   Вычисление невязок
  allocate(y1(n));
  y1 = matmul(In,x);
  E = y - y1;
  !  E=y-matmul(In,x) ! почему то при вычислении этой строчки в лоб тратиться очень много памяти и мы имеет StackOverflow. Поэтому сделано через динамический массив
  E1=dot_product(E,E);
  deallocate(y1);
  
  !---- Вычисление среднеквадратического отклонения по невязкам (критерий хи-квадрат)
  s0_dis=sqrt(E1/(n-m))
  
  !---- вычисление ошибок коэффициентов решения
  do i=1,m
  sigma(i)= sqrt(Cov(i,i));
  enddo 
  
  !---- Вычисление корреляционной матрицы
  do i=1,m
  do j=1,m
  Corr(i,j)=Cov(i,j)/sqrt(Cov(i,i)*Cov(j,j))
  enddo
  enddo;
  
  endsubroutine TLS_Gen;
  
}