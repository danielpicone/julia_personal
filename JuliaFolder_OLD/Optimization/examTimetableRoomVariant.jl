# Exam Timetable optimisation
# This is a script to minimise the number of days in an exam timetable
# Each student i has 4 exams, with the objective to minimise the number of days
# of exams

using JuMP
using CPLEX
start = time()
timetable = Model(solver = CplexSolver())


# Let i = {1,2,...N} be the students
N = 10
# Let j = {1,2,...,M} be the exams
M = 5
# Let k = {1,2,...,K} be the slots for exams to occur in
K = 100
# Let l = {1,...,L} be the number of rooms available per day
L = 10
# Maximum number of exams a student may have
max_exams = 4
r_min = 5

@variable(timetable, x[j=1:M, k=1:K, l=1:L],Bin)
@variable(timetable, r[j=1:M, l=1:L], Bin)
@variable(timetable, y[k=1:K], Bin)
# Create the constraint that ensures all exams occur once only
for j=1:M
    @constraint(timetable, sum(x[j,k,l] for k=1:K,l=1:L) == 1)
end
# Create the constraint that slot k-1 is filled before slot k
for k=2:K
    @constraint(timetable, y[k] <= y[k-1])
end
# Create the constraint that ensures a slot is active if there is an exam in it
for k=1:K
    @constraint(timetable, y[k] >= sum(x[j,k,l] for j=1:M,l=1:L))
end
# Create the constraint that every student must have their exams in a different
# slot
# First randomly assign student their exams
student = Array{Int64}(N,max_exams)
# Want a distribution of:
# 80% have 4 exams
E1 = 1.0
# 10% have 3 exams
E2 = E1 + 0.1
# 5% have 2 exams
E3 = E2 + 0.05
# 5% have 1 exam
E4 = E3 + 0.05

for i=1:N
    if rand() <= E1
        student[i,:] = randperm(M)[1:max_exams]
    elseif E1 < rand() <= E2
        student[i,:] = [randperm(M)[1:max_exams-1]; K+1]
    elseif E2 < rand() <= E3
        student[i,:] = [randperm(M)[1:max_exams-2]; K+1; K+1]
    else
        student[i,:] = [randperm(M)[1]; K+1; K+1; K+1]
    end
end

max_exam_size = 0
for j=1:M
    l = 1
    max_exam_size = max(count(l->l==j,student[:,:]),max_exam_size)
end

# Create the room sizes, which go from r_min to max_exam_size
c = Array{Int64}(L)
for l=1:L
    c[l] = rand(r_min:max_exam_size+10)
end
# Create the constraint which ensures that every exam is in a room big enough
for j=1:M
    for l=1:L
        @constraint(timetable, sum(x[j,k,l] for k = 1:K) <= c[l]*r[j,l])
    end
end

# Now create the constraint ensuring every exam sat by student i must be on a different day
for k=1:K
    for i=1:N
        p = 1
        count(p->p<=K,student[i,:])
        @constraint(timetable, sum(x[j,k,l] for l=1:L,j in student[i,1:p])<= 1)
    end
end

# Add the objective function
@objective(timetable, Min, sum(y[k] for k=1:K))

solve(timetable)

exam_slots = getvalue(y)
num_slots = ceil(sum(exam_slots[k] for k=1:K)/2)
@printf("The number of days = %d\n",num_slots)

# To see when student i has exams
getvalue(r)
# Test if all students in exam j can fit in the room assigned
for j=1:M
    # Find the room the exam is in
    #find(p->p==1,getvalue(r)[j,:]) <=
end


"""stud = 10
#l = 1
l1 = count(l1->l1<=K,student[stud,:])
days = Array{Int64}(l1)
k = 1
for j in student[stud,1:l1]
    days[k] = ceil(find(x1->x1>=0.9999999,getvalue(x[j,:]))[1]/2)
    k += 1
end
@printf("Student %d has exams on days: ",stud)
for i=1:length(days)
    println(days[i])
end"""
