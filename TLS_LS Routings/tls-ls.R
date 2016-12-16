TLS_Gen <- function (mode, scaling, ef, n, m, In, y, X, E, S, Corr, sigma, s0, s0_dis, err)
{
  # m - variables qty
  # n - equations qty
  # mode - 
  # scaling - sort of scaling in case of TLS solution, bit mask: 1 - by rows, 2 - by columns, 3 - Frobenian
  # ef - EIV qty
  # In - equations
  # y - observations
  # X - solution
  # E - discripency 
  # S - Singular Values
  # Corr - covariance matrix
  # sigma - solution errors
  # s0 - mean square deviation
  # 
  # err - error code
  
  
  #---------------------------------------------------------------
  #     (TLS, МНК) в общем виде
  
  #integer, intent(in) :: mode, scaling, ef, n, m
  #doubleprecision, intent(in) :: In(n,m),y(n),p(n)
  #doubleprecision, intent(out):: x(m), E(n), S(m+1), Corr(m,m), sigma(m), s0, s0_dis, err;  
  
  #doubleprecision, allocatable :: a(:,:), b(:,:), y1(:), U(:,:), V(:,:), dn(:,:), cn(:,:), 
  #doubleprecision, allocatable :: b11(:,:), b12(:,:), b21(:,:), b22(:,:), b_shur(:,:), b_inv(:,:)  
  #doubleprecision, allocatable :: b12v(:), b21v(:)
  #integer i,j, nef
  #doubleprecision ee(m,m),c(m,m),d(m,m),d1(m,m),det,E1, b11s, Cov(m,m);

  
  #complex(kind(1d0)) s_c(m+1) 
  #complex(kind(1d0)), allocatable :: s_b(:)
  
  if ((m>n) | (m==n) | (n<3) | (m<2))
  { 
    err<<- (-65535)
    return;
  }
  
  
  S<-vector("numeric", m+1) ; X<-vector("numeric", m); E<-vector("numeric", n); Corr<-matrix(0, m, m); sigma<-vector("numeric", m); s0<-0; s0_dis<-0;
  
  # ===================================================================
  # ========= вычисление сингулярных чисел матрицы исходных данных для оценки числа обусловленности  
  # ========= The Total Least Square Problem, p. 232
  # ===================================================================
  in_svd <- svd(In)  

  # =========== Подготовка исходных данных, вычисляем дополненную матрицу условных уравнений
  
  a <- cbind(In, y);
  
  # ===================================================================
  # ========= Нормализация исходных данных
  # ========= The Total Least Square Problem, p. 90, $3.6.2, (3.125)
  # ===================================================================
  if (scaling>0)
  {
  
      dn <- matrix(0, n, n);
      cn <- matrix(0, m+1, m+1);
      
      # по строчкам   
      
      if ((scaling==1) || (scaling==3))
      {
        for (i in 1:n)
          dn[i,i] <- 1 / sqrt(sum((a[i,]^2)));
                              
        a <- dn %*% a;
      }
      
      # по столбцам
      
      if ((scaling==2) || (scaling==3)) 
      {
         for(i in 1:m+1)
           cn[i,i] <- 1 / sqrt(sum((a[,i]^2)));
         
        a <- a %*% cn
      }
      
      # норма Фробениуса
      #if (scaling==4){ 
      #  call NR2RR (a, fn)   #  !!!
      #}
      
  }
  
  # ===================================================================
  # ====================   Вычисление решения
  # ============ Решение через сингуляргое разложение
  # ===================================================================
  if (mode == 1){
     a_svd <- svd(a); 
    
     #----  проверка решения
     S <- a_svd$d  # вектор сингулярных чисел
     aa <- a_svd$u %*% (diag(a_svd$d) %*% t(a_svd$v) );  
     ee <- aa-a; 
     
     #----  вычисляем х
     x2 <-  a_svd$v[,m+1];  
     
     if ((scaling==2) || (scaling==3)) 
        x2 = cn %*% x2; 
     
     x2 = -1*x2/x2[m+1]
     X = x2[1:m];
     
     #---- Вычисление числа обусловленности
     #---- вот это надо проверить, Branham, Astronomical Data Reduction with TLS, p.655
     err<- (in_svd$d[1]^2-a_svd$d[m+1]^2)/(in_svd$d[m]^2-a_svd$d[m+1]^2);
     
     #----  вычисление дисперсию, The Total Least Squares Problem, 8.15
     s0 <- a_svd$d[m+1] / sqrt(n+0.0);
     
     #---- вычисление ковариационной матрицы, The Total Least Squares Problem, page 242, 8.47
     d <- matrix(0, m, m) ;
     diag(d) <- n*a_svd$d[m+1]^2;
     d1 <- (t(In) %*% In) - d;
     Cov <- solve(d1)
     Cov <- (1+sum(X^2)) * s0^2 * Cov;
     
  } else if (mode == 2) {
  # ===================================================================
  # ===========  Решение через собственные числа
  # ===================================================================
      ss <- t(a) %*% a;  # номальная система уравнений
  
      if ((ef<1) || (ef>=m))# вариант польного TLS, все переменные содержат ошибки
      {
          nef <- m+1;
      
          s_c <- eigen(ss);
          S <- s_c$val;
          
          if (S[m+1]<0) S[m+1] <- 0;
          ee <- diag(rep(S[m+1], m))

          #----  вычисление sigma_0, The Total Least Squares Problem, 8.15
          s0 <- sqrt(S[m+1])/sqrt(n+0.0);
          #---- Вычисление числа обусловленности
          err<- (in_svd$d[1]^2-S[m+1])/(in_svd$d[m]^2-S[m+1]);
          
          #---- неполное вычисление ковариационной матрицы, по аналогии с The Total Least Squares Problem, page 242, 8.47
          
          d <- matrix(0, m , m);
          d <- diag(rep(S[m+1], m));
          d1 <- t(In) %*% In - d;
          Cov <- solve(d1);
          
          
  #    этот частный случай (когда только первый столбец error-free) вошел в более общий, реализованный ниже. Код рабочий, на всякий случай оставляю
  #    elseif (ef==-1) then     
  #     nef=m;
  #     allocate (b21v(m),b12v(m),b22(m, m), b_shur(m, m), s_b(m));
  #     b21v = ss(1,2:m+1);
  #     b12v = ss(2:m+1, 1);
  #     b22 = ss(2:m+1,2:m+1);
  #     do i=1,m
  #       do j=1,m
  #         b_shur(i,j) = b21v(i)*b12v(j)
  #       enddo
  #     enddo
  
  #     b_shur = b_shur/ss(1,1);
  #     b_shur = b22 - b_shur
  
  #!     call lin_eig_gen(b_shur, s_b);     
  #!     ee = 0;
  #!     s=0;
  #!     do i=2,m
  #!       ee(i,i) = s_b(m);      
  #!       s(i-1) = s_b(i-1)
  #!     enddo   
  #!     s(m) = s_b(m);
  #!     deallocate(b21v, b12v, b22, b_shur, s_b);
  
      } else {   #! вариант TLS-LS, когда первый ef-столбцов не содержат ошибки, а последние nef - содержат
          nef <- m+1-ef; 
          
          #allocate (b21(nef,ef), b12(ef,nef),b22(nef, nef), b_shur(nef, nef),b_inv(ef,ef), s_b(nef), b11(ef,ef));
          
          b12 <- as.matrix(ss[1:ef,(ef+1):(m+1)]);
          if (ef == 1) b12 <- t(b12)
          b21 <- as.matrix(ss[(ef+1):(m+1), 1:ef]);
          if (ef == m) b21 <- t(b21)
          b22 <- as.matrix(ss[(ef+1):(m+1),(ef+1):(m+1)]);
          b11 <- as.matrix(ss[1:ef,1:ef]);
          
          b_inv <- solve(b11)

          b_shur <- b22 - (b21 %*% b_inv) %*% b12;
          
          s_b <- eigen(b_shur)
          
          #ee <- matrix(0, m, m)
          #for (i in ef+1:m)
          #  ee(i,i) <- s_b$val[nef];
          
          ee <- diag(c(rep(0, ef), rep(s_b$val[nef], m-ef)))
          
          #S <- vector("numeric", m+1);
          #for (i in 1:nef)
          #  S[i] = s_b$val[i]
          S[1:nef] <- s_b$val
          
          #----  вычисление sigma_0, по аналогии с The Total Least Squares Problem, 8.15, возможно не правильно, надо как-то проверить
          s0 <- sqrt(S[nef])/sqrt(n+0.0);
          
          #---- неполное вычисление ковариационной матрицы, по аналогии с The Total Least Squares Problem, page 242, 8.47, возможно не правильно, надо как-то проверить
          
          d <- matrix(0, m , m);
          d <- diag(rep(s_b$val[nef], m));
          
          d1 <- t(In) %*% In - d;
          Cov <- solve(d1);
          
          #---- Вычисление числа обусловленности
          err<- (in_svd$d[1]^2-S[nef])/(in_svd$d[m]^2-S[nef]);
  
      }
  
    #--------- Вычисление матрицы нормальных уравнений
    c<-t(a[,1:m]) %*% a[,1:m]-ee;  # случай с масштабированием
    
    #---------  вычисление обратной матрицы
    d <- solve(c) 
    ee <- c %*% d;  # проверка
    
    #----   вычисление решения  
    X <- as.vector(d %*% (t(a[,1:m]) %*% a[,m+1])); # с учетом масштабирования
    
    if ((scaling==2) || (scaling==3))
        X = cn %*% x / cn(m+1,m+1); 
    
    #---- заканчиваем вычисление корреляции
    Cov = (1+sum(X^2)) * s0^2 * Cov;
    
  }


  #-------   Вычисление невязок
  E <- y - In %*% X;
  
  #---- Вычисление среднеквадратического отклонения по невязкам 
  s0_dis=sqrt(sum(E*E)/(n-m))
  
  #---- вычисление ошибок коэффициентов решения
  sigma <- sqrt(diag(Cov))
  
  
  #---- Вычисление корреляционной матрицы
  Corr <- matrix(0, m, m)
  for (i in 1:m)
    for (j in 1:m)
      Corr[i,j]=Cov[i,j]/sqrt(Cov[i,i]*Cov[j,j])
  
  print(X);
  print(s0);
  print(s0_dis);
  print(sigma);
  print(err)
  #print(Cov);

}