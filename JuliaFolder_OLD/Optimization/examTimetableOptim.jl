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
# We will create simple_students and then duplicate them and change 1 exam each
repeated_students = 10
simple_students = Int(N/repeated_students)
# Let j = {1,2,...,M} be the exams
M = 5
# Let k = {1,2,...,K} be the slots for exams to occur in
K = 5
# Maximum number of exams a student may have
max_exams = 4

@variable(timetable, x[j=1:M,k=1:K],Bin)
@variable(timetable, y[k=1:K], Bin)
# Create the constraint that ensures all exams occur once only
for j=1:M
    @constraint(timetable, sum(x[j,k] for k=1:K) == 1)
end
# Create the constraint that slot k-1 is filled before slot k
for k=2:K
    @constraint(timetable, y[k] <= y[k-1])
end
# Create the constraint that ensures a slot is active if there is an exam in it
for k=1:K
    @constraint(timetable, y[k] >= sum(x[j,k] for j=1:M))
end
# Create the constraint that every student must have their exams in a different
# slot
# First randomly assign student their exams
student = Array{Int64}(N,max_exams)
# Want a distribution of:
# 80% have 4 exams
E1 = 0.8
# 10% have 3 exams
E2 = E1 + 0.1
# 5% have 2 exams
E3 = E2 + 0.05
# 5% have 1 exam
E4 = E3 + 0.05

for i=1:simple_students
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

full_student = Array{Int64}(N,4)
for i=1:simple_students
    for j=1:repeated_students
        student[i,1] = randperm(M)[1]
        full_student[j+(i-1)*repeated_students,:] = student[i,:]
    end
end

student = full_student

# Now create the constraint
for k=1:K
    for i=1:N
        l = 1
        count(l->l<=K,student[i,:])
        @constraint(timetable, sum(x[j,k] for j in student[i,1:l])<= 1)
    end
end

# Add the objective function
@objective(timetable, Min, sum(y[k] for k=1:K))

solve(timetable)

exam_slots = getvalue(y)
num_slots = ceil(sum(exam_slots[k] for k=1:K)/2)
@printf("The number of days = %d\n",num_slots)

# To see when student i has exams
stud = 10
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
end
