# Creating role
aws iam create-role --role-name lambda-vpc-role --assume-role-policy-document file://trust-policy.json
aws iam attach-role-policy --role-name lambda-vpc-role --policy-arn arn:aws:iam::aws:policy/service-role/AWSLambdaBasicExecutionRole
aws iam attach-role-policy --role-name lambda-vpc-role --policy-arn arn:aws:iam::aws:policy/service-role/AWSLambdaVPCAccessExecutionRole

# Creating EC2 VM (group -> default VPC security group, incoming UDP connections on port 10,000 allowed)
aws ec2 run-instances --image-id ami-06e4ca05d431835e9 --count 1 --instance-type t2.micro --key-name MyKeyPair --security-group-ids sg-06652ad96a1727a6a --subnet-id subnet-0a105180830223542 
#
## To get the IP
#aws ec2 describe-instances

# Creating function
aws lambda create-function \
--function-name  \
--role arn:aws:iam::183425780977:role/lambda-vpc-role \
--runtime provided \
--timeout 2 \
--memory-size 512 \
--handler fmi_svd_bench \
--zip-file fileb://fmi_svd_bench.zip

# Compiling and packaging code into zip file
make aws-lambda-package-fmi_svd_bench

# Updating code with newest zip
#aws lambda update-function-code --function-name smibenchmark --zip-file fileb://smibenchmark.zip
